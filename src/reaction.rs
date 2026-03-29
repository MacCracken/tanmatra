//! Nuclear reactions and Q-value calculations.
//!
//! Provides nuclear reaction types, Q-value computation from mass differences,
//! Coulomb barrier estimates, and preset reactions for common fusion and fission
//! processes.

use crate::constants::COULOMB_MEV_FM;
use crate::error::TanmatraError;
use crate::nucleus::Nucleus;
use alloc::string::String;
use serde::{Deserialize, Serialize};

/// A nuclear reaction: reactants -> products + Q.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NuclearReaction {
    /// Name or description of the reaction.
    pub name: String,
    /// Projectile nucleus.
    pub projectile: Nucleus,
    /// Target nucleus.
    pub target: Nucleus,
    /// Product nuclei.
    pub products: alloc::vec::Vec<Nucleus>,
    /// Q-value in MeV (positive = exothermic, negative = endothermic).
    pub q_value_mev: f64,
}

/// Calculates the Q-value of a reaction from the mass difference.
///
/// Q = (sum of reactant masses - sum of product masses) * c^2
///
/// All masses in MeV/c^2. Positive Q means energy is released.
#[must_use]
#[inline]
pub fn q_value(reactant_masses: &[f64], product_masses: &[f64]) -> f64 {
    let sum_reactants: f64 = reactant_masses.iter().sum();
    let sum_products: f64 = product_masses.iter().sum();
    sum_reactants - sum_products
}

/// Estimates the Coulomb barrier for two nuclei approaching each other.
///
/// V_c = k_e * Z1 * Z2 * e^2 / (R1 + R2)
///
/// where R = r0 * A^(1/3) and k_e*e^2 = 1.44 MeV*fm.
///
/// Returns the barrier height in MeV.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidReaction`] if either nucleus is invalid.
#[inline]
pub fn coulomb_barrier(n1: &Nucleus, n2: &Nucleus) -> Result<f64, TanmatraError> {
    let r1 = n1.nuclear_radius();
    let r2 = n2.nuclear_radius();
    let r_sum = r1 + r2;

    if r_sum <= 0.0 {
        return Err(TanmatraError::InvalidReaction(String::from(
            "nuclear radii sum to zero",
        )));
    }

    let z1 = n1.z() as f64;
    let z2 = n2.z() as f64;

    Ok(COULOMB_MEV_FM * z1 * z2 / r_sum)
}

// ---------------------------------------------------------------------------
// Cross-section calculations
// ---------------------------------------------------------------------------

/// Returns the geometric cross-section of a nucleus in barns (1 barn = 1e-24 cm^2).
///
/// σ_geo = π R^2 where R = r0 * A^(1/3) in femtometers.
///
/// 1 barn = 100 fm^2.
#[must_use]
#[inline]
pub fn geometric_cross_section_barns(nucleus: &Nucleus) -> f64 {
    let r = nucleus.nuclear_radius(); // in fm
    let sigma_fm2 = core::f64::consts::PI * r * r;
    sigma_fm2 / 100.0 // convert fm^2 to barns
}

/// Breit-Wigner single-level resonance cross-section.
///
/// σ(E) = π λ̄² g (Γ_i Γ_f) / ((E - E_r)² + (Γ/2)²)
///
/// where:
/// - `energy_mev`: projectile kinetic energy in MeV
/// - `resonance_energy_mev`: resonance energy E_r in MeV
/// - `total_width_mev`: total decay width Γ in MeV
/// - `partial_in_mev`: entrance channel partial width Γ_i
/// - `partial_out_mev`: exit channel partial width Γ_f
/// - `g_factor`: statistical spin factor g = (2J+1)/((2s1+1)(2s2+1))
/// - `reduced_mass_mev`: reduced mass of projectile-target system in MeV/c²
///
/// Returns cross-section in barns.
#[must_use]
#[inline]
pub fn breit_wigner_cross_section(
    energy_mev: f64,
    resonance_energy_mev: f64,
    total_width_mev: f64,
    partial_in_mev: f64,
    partial_out_mev: f64,
    g_factor: f64,
    reduced_mass_mev: f64,
) -> f64 {
    if energy_mev <= 0.0 || reduced_mass_mev <= 0.0 {
        return 0.0;
    }

    // de Broglie wavelength squared: λ̄² = ħ²c² / (2 μ E)
    // ħc = 197.3269804 MeV·fm, so ħ²c² = 197.327² MeV²·fm²
    let hbar_c = 197.326_980_4; // MeV·fm
    let lambda_bar_sq = hbar_c * hbar_c / (2.0 * reduced_mass_mev * energy_mev);

    let de = energy_mev - resonance_energy_mev;
    let half_width = total_width_mev / 2.0;
    let denominator = de * de + half_width * half_width;

    if denominator <= 0.0 {
        return 0.0;
    }

    let sigma_fm2 =
        core::f64::consts::PI * lambda_bar_sq * g_factor * partial_in_mev * partial_out_mev
            / denominator;

    // Convert fm² to barns (1 barn = 100 fm²)
    sigma_fm2 / 100.0
}

/// Known thermal neutron cross-sections (at 0.0253 eV) in barns.
///
/// Source: NNDC, Mughabghab "Atlas of Neutron Resonances" (2018).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ThermalCrossSection {
    /// Target nucleus.
    pub target: Nucleus,
    /// Thermal neutron absorption cross-section in barns.
    pub absorption_barns: f64,
    /// Thermal neutron fission cross-section in barns (0 if not fissile).
    pub fission_barns: f64,
    /// Thermal neutron scattering cross-section in barns.
    pub scattering_barns: f64,
}

/// Returns thermal neutron cross-sections for key isotopes.
///
/// Source: NNDC evaluated nuclear data, Mughabghab (2018).
#[must_use]
pub fn thermal_neutron_cross_sections() -> alloc::vec::Vec<ThermalCrossSection> {
    alloc::vec![
        // U-235: fissile
        ThermalCrossSection {
            target: Nucleus::uranium_235(),
            absorption_barns: 680.9,
            fission_barns: 585.1,
            scattering_barns: 15.04,
        },
        // U-238: fertile
        ThermalCrossSection {
            target: Nucleus::uranium_238(),
            absorption_barns: 2.680,
            fission_barns: 0.0,
            scattering_barns: 8.871,
        },
        // Pu-239: fissile
        ThermalCrossSection {
            target: Nucleus::new(94, 239).unwrap_or_else(|_| Nucleus::hydrogen_1()),
            absorption_barns: 1011.3,
            fission_barns: 747.4,
            scattering_barns: 7.7,
        },
        // B-10: neutron absorber
        ThermalCrossSection {
            target: Nucleus::new(5, 10).unwrap_or_else(|_| Nucleus::hydrogen_1()),
            absorption_barns: 3835.0,
            fission_barns: 0.0,
            scattering_barns: 2.16,
        },
        // Cd-113: strong absorber
        ThermalCrossSection {
            target: Nucleus::new(48, 113).unwrap_or_else(|_| Nucleus::hydrogen_1()),
            absorption_barns: 20_600.0,
            fission_barns: 0.0,
            scattering_barns: 12.1,
        },
        // H-1: moderator
        ThermalCrossSection {
            target: Nucleus::hydrogen_1(),
            absorption_barns: 0.3326,
            fission_barns: 0.0,
            scattering_barns: 82.02,
        },
        // C-12: moderator
        ThermalCrossSection {
            target: Nucleus::carbon_12(),
            absorption_barns: 0.003_53,
            fission_barns: 0.0,
            scattering_barns: 5.551,
        },
        // Fe-56: structural material
        ThermalCrossSection {
            target: Nucleus::iron_56(),
            absorption_barns: 2.59,
            fission_barns: 0.0,
            scattering_barns: 12.42,
        },
    ]
}

// --- Preset reactions with real Q-values ---

/// D-T fusion: D + T -> He-4 + n, Q = 17.6 MeV.
///
/// The most favorable fusion reaction for energy production.
#[must_use]
pub fn dt_fusion() -> NuclearReaction {
    NuclearReaction {
        name: String::from("D-T fusion"),
        projectile: Nucleus::new(1, 2).unwrap_or_else(|_| Nucleus::hydrogen_1()),
        target: Nucleus::new(1, 3).unwrap_or_else(|_| Nucleus::hydrogen_1()),
        products: alloc::vec![Nucleus::helium_4()],
        // He-4 + neutron; Q = 17.6 MeV
        q_value_mev: 17.6,
    }
}

/// D-D fusion (He-3 branch): D + D -> He-3 + n, Q = 3.27 MeV.
#[must_use]
pub fn dd_fusion_he3() -> NuclearReaction {
    NuclearReaction {
        name: String::from("D-D fusion (He-3 branch)"),
        projectile: Nucleus::new(1, 2).unwrap_or_else(|_| Nucleus::hydrogen_1()),
        target: Nucleus::new(1, 2).unwrap_or_else(|_| Nucleus::hydrogen_1()),
        products: alloc::vec![Nucleus::new(2, 3).unwrap_or_else(|_| Nucleus::helium_4())],
        q_value_mev: 3.27,
    }
}

/// D-D fusion (tritium branch): D + D -> T + p, Q = 4.03 MeV.
#[must_use]
pub fn dd_fusion_t() -> NuclearReaction {
    NuclearReaction {
        name: String::from("D-D fusion (tritium branch)"),
        projectile: Nucleus::new(1, 2).unwrap_or_else(|_| Nucleus::hydrogen_1()),
        target: Nucleus::new(1, 2).unwrap_or_else(|_| Nucleus::hydrogen_1()),
        products: alloc::vec![
            Nucleus::new(1, 3).unwrap_or_else(|_| Nucleus::hydrogen_1()),
            Nucleus::hydrogen_1(),
        ],
        q_value_mev: 4.03,
    }
}

/// Proton-proton chain step 1: p + p -> D + e+ + neutrino, Q = 0.42 MeV.
///
/// The dominant energy source in the Sun.
#[must_use]
pub fn pp_chain_step1() -> NuclearReaction {
    NuclearReaction {
        name: String::from("p-p chain (step 1)"),
        projectile: Nucleus::hydrogen_1(),
        target: Nucleus::hydrogen_1(),
        products: alloc::vec![Nucleus::new(1, 2).unwrap_or_else(|_| Nucleus::hydrogen_1())],
        // D + positron + neutrino; Q ≈ 0.42 MeV
        q_value_mev: 0.42,
    }
}

/// Uranium-235 fission (typical): U-235 + n -> Ba-141 + Kr-92 + 3n, Q ≈ 200 MeV.
///
/// One of many possible fission channels. The average energy release per
/// fission of U-235 is approximately 200 MeV.
#[must_use]
pub fn u235_fission() -> NuclearReaction {
    NuclearReaction {
        name: String::from("U-235 fission (typical)"),
        projectile: Nucleus::new(1, 1).unwrap_or_else(|_| Nucleus::hydrogen_1()), // neutron approximated as H-1
        target: Nucleus::uranium_235(),
        products: alloc::vec![
            Nucleus::new(56, 141).unwrap_or_else(|_| Nucleus::iron_56()), // Ba-141
            Nucleus::new(36, 92).unwrap_or_else(|_| Nucleus::iron_56()),  // Kr-92
        ],
        // + 3 neutrons; Q ≈ 200 MeV
        q_value_mev: 200.0,
    }
}

/// CNO cycle (net): 4p -> He-4 + 2e+ + 2nu, Q = 25.03 MeV.
///
/// The dominant energy source in stars more massive than ~1.3 solar masses.
#[must_use]
pub fn cno_cycle() -> NuclearReaction {
    NuclearReaction {
        name: String::from("CNO cycle (net)"),
        projectile: Nucleus::hydrogen_1(),
        target: Nucleus::carbon_12(), // Catalyst
        products: alloc::vec![Nucleus::helium_4(), Nucleus::carbon_12()],
        q_value_mev: 25.03,
    }
}

/// Triple-alpha process: 3 He-4 -> C-12, Q = 7.275 MeV.
///
/// The process by which stars synthesize carbon from helium.
#[must_use]
pub fn triple_alpha() -> NuclearReaction {
    NuclearReaction {
        name: String::from("Triple-alpha process"),
        projectile: Nucleus::helium_4(),
        target: Nucleus::helium_4(),
        products: alloc::vec![Nucleus::carbon_12()],
        q_value_mev: 7.275,
    }
}

// ---------------------------------------------------------------------------
// Nucleosynthesis pathways
// ---------------------------------------------------------------------------

/// A nucleosynthesis process type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum NucleosynthesisProcess {
    /// Slow neutron capture (s-process): captures slower than beta-decay.
    SProcess,
    /// Rapid neutron capture (r-process): captures faster than beta-decay.
    RProcess,
    /// Proton capture (rp-process).
    RpProcess,
}

/// A step in a nucleosynthesis pathway.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NucleosynthesisStep {
    /// The nucleus at this step.
    pub nucleus: Nucleus,
    /// The reaction that produced this nucleus (neutron capture, beta decay, etc.).
    pub reaction: String,
}

/// A nucleosynthesis pathway: a sequence of nuclear reactions building heavier elements.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NucleosynthesisPathway {
    /// Name of the pathway.
    pub name: String,
    /// The process type.
    pub process: NucleosynthesisProcess,
    /// Seed nucleus.
    pub seed: Nucleus,
    /// Steps in the pathway.
    pub steps: alloc::vec::Vec<NucleosynthesisStep>,
}

/// Returns the main s-process pathway from Fe-56 seed.
///
/// The s-process (slow neutron capture) occurs in AGB stars. Neutrons are
/// captured one at a time, with beta-decay occurring between captures when
/// the nucleus becomes unstable. This traces the valley of stability.
///
/// The pathway shown is a simplified main component through key s-process
/// abundance peaks at Sr-88, Ba-138, and Pb-208.
///
/// Source: B2FH (1957), Kaeppeler et al. Rev. Mod. Phys. 83, 157 (2011).
#[must_use]
pub fn s_process_main() -> NucleosynthesisPathway {
    NucleosynthesisPathway {
        name: String::from("s-process main component"),
        process: NucleosynthesisProcess::SProcess,
        seed: Nucleus::iron_56(),
        steps: alloc::vec![
            NucleosynthesisStep {
                nucleus: Nucleus::new(26, 57).unwrap_or_else(|_| Nucleus::iron_56()),
                reaction: String::from("Fe-56(n,γ)Fe-57"),
            },
            NucleosynthesisStep {
                nucleus: Nucleus::new(26, 58).unwrap_or_else(|_| Nucleus::iron_56()),
                reaction: String::from("Fe-57(n,γ)Fe-58"),
            },
            NucleosynthesisStep {
                nucleus: Nucleus::new(27, 59).unwrap_or_else(|_| Nucleus::iron_56()),
                reaction: String::from("Fe-59→Co-59 (β⁻)"),
            },
            // Skip ahead to key s-process peaks
            NucleosynthesisStep {
                nucleus: Nucleus::new(38, 88).unwrap_or_else(|_| Nucleus::iron_56()),
                reaction: String::from("...→Sr-88 (1st peak, N=50)"),
            },
            NucleosynthesisStep {
                nucleus: Nucleus::new(56, 138).unwrap_or_else(|_| Nucleus::iron_56()),
                reaction: String::from("...→Ba-138 (2nd peak, N=82)"),
            },
            NucleosynthesisStep {
                nucleus: Nucleus::new(82, 208).unwrap_or_else(|_| Nucleus::iron_56()),
                reaction: String::from("...→Pb-208 (3rd peak, N=126)"),
            },
        ],
    }
}

/// Returns the main r-process pathway.
///
/// The r-process (rapid neutron capture) occurs in neutron star mergers
/// and core-collapse supernovae. Neutrons are captured much faster than
/// beta-decay, driving nuclei far from stability toward the neutron drip
/// line before beta-decaying back to stability.
///
/// Key r-process abundance peaks occur at A ≈ 80, 130, and 195, which
/// correspond to neutron magic numbers N = 50, 82, 126 on the neutron-rich
/// side of stability.
///
/// Source: Cowan et al., Rev. Mod. Phys. 93, 015002 (2021).
#[must_use]
pub fn r_process_main() -> NucleosynthesisPathway {
    NucleosynthesisPathway {
        name: String::from("r-process main component"),
        process: NucleosynthesisProcess::RProcess,
        seed: Nucleus::iron_56(),
        steps: alloc::vec![
            // r-process peaks (after beta-decay back to stability)
            NucleosynthesisStep {
                nucleus: Nucleus::new(34, 80).unwrap_or_else(|_| Nucleus::iron_56()),
                reaction: String::from("1st peak: Se-80 (A≈80, from N=50 waiting point)"),
            },
            NucleosynthesisStep {
                nucleus: Nucleus::new(52, 130).unwrap_or_else(|_| Nucleus::iron_56()),
                reaction: String::from("2nd peak: Te-130 (A≈130, from N=82 waiting point)"),
            },
            NucleosynthesisStep {
                nucleus: Nucleus::new(76, 195).unwrap_or_else(|_| Nucleus::iron_56()),
                reaction: String::from("3rd peak: Os/Pt-195 (A≈195, from N=126 waiting point)"),
            },
            NucleosynthesisStep {
                nucleus: Nucleus::new(92, 238).unwrap_or_else(|_| Nucleus::uranium_238()),
                reaction: String::from("Actinide production: U-238"),
            },
        ],
    }
}

// ---------------------------------------------------------------------------
// Neutron moderation and thermalization
// ---------------------------------------------------------------------------

/// Calculates the maximum fractional energy loss per elastic collision.
///
/// For a neutron (mass 1) scattering off a nucleus of mass number A:
/// α = ((A-1)/(A+1))², and the maximum fractional loss = 1 - α.
///
/// Returns the fraction of energy lost in a single head-on collision.
#[must_use]
#[inline]
pub fn max_energy_loss_fraction(target_a: u32) -> f64 {
    let a = target_a as f64;
    let alpha = ((a - 1.0) / (a + 1.0)) * ((a - 1.0) / (a + 1.0));
    1.0 - alpha
}

/// Calculates the average logarithmic energy decrement per collision (ξ).
///
/// ξ = 1 + ((A-1)²/(2A)) * ln((A-1)/(A+1))
///
/// For hydrogen (A=1), ξ = 1.0 exactly.
/// This is the average lethargy gain per collision.
#[must_use]
#[inline]
pub fn average_lethargy_gain(target_a: u32) -> f64 {
    if target_a <= 1 {
        return 1.0; // hydrogen: ξ = 1
    }
    let a = target_a as f64;
    // ξ = 1 + ((A-1)²/(2A)) * ln((A-1)/(A+1))
    let am1 = a - 1.0;
    let ap1 = a + 1.0;
    1.0 + (am1 * am1 / (2.0 * a)) * libm::log(am1 / ap1)
}

/// Calculates the number of collisions needed to thermalize a neutron.
///
/// n = ln(E_initial / E_final) / ξ
///
/// Typical: E_initial = 2 MeV (fission neutron), E_final = 0.0253 eV (thermal).
///
/// Returns the number of collisions (as f64 since it's an average).
#[must_use]
#[inline]
pub fn collisions_to_thermalize(target_a: u32, e_initial_ev: f64, e_final_ev: f64) -> f64 {
    let xi = average_lethargy_gain(target_a);
    if xi <= 0.0 || e_initial_ev <= 0.0 || e_final_ev <= 0.0 {
        return 0.0;
    }
    libm::log(e_initial_ev / e_final_ev) / xi
}

/// Calculates the moderating ratio (ξ Σ_s / Σ_a).
///
/// A higher moderating ratio means better moderator performance.
/// Typical values: H₂O ≈ 72, D₂O ≈ 5670, graphite ≈ 192.
///
/// Parameters:
/// - `lethargy_gain`: ξ value
/// - `scatter_xs`: macroscopic scattering cross-section (cm⁻¹)
/// - `absorb_xs`: macroscopic absorption cross-section (cm⁻¹)
#[must_use]
#[inline]
pub fn moderating_ratio(lethargy_gain: f64, scatter_xs: f64, absorb_xs: f64) -> f64 {
    if absorb_xs <= 0.0 {
        return f64::INFINITY;
    }
    lethargy_gain * scatter_xs / absorb_xs
}

// ---------------------------------------------------------------------------
// Fission product yields
// ---------------------------------------------------------------------------

/// A fission product yield data point.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct FissionYield {
    /// Mass number A of the fission product.
    pub mass_number: u32,
    /// Cumulative fission yield (fraction per fission).
    pub yield_fraction: f64,
}

/// Returns thermal neutron fission yield distribution for U-235.
///
/// These are cumulative yields at key mass numbers showing the
/// characteristic double-humped distribution.
///
/// Source: ENDF/B-VIII.0 fission yield data.
#[must_use]
pub fn u235_fission_yields() -> alloc::vec::Vec<FissionYield> {
    alloc::vec![
        FissionYield {
            mass_number: 85,
            yield_fraction: 0.0131
        },
        FissionYield {
            mass_number: 90,
            yield_fraction: 0.0580
        },
        FissionYield {
            mass_number: 95,
            yield_fraction: 0.0650
        },
        FissionYield {
            mass_number: 99,
            yield_fraction: 0.0611
        },
        FissionYield {
            mass_number: 101,
            yield_fraction: 0.0514
        },
        FissionYield {
            mass_number: 105,
            yield_fraction: 0.0093
        },
        FissionYield {
            mass_number: 110,
            yield_fraction: 0.0002
        },
        FissionYield {
            mass_number: 115,
            yield_fraction: 0.0001
        },
        FissionYield {
            mass_number: 120,
            yield_fraction: 0.0001
        },
        FissionYield {
            mass_number: 131,
            yield_fraction: 0.0290
        },
        FissionYield {
            mass_number: 133,
            yield_fraction: 0.0670
        },
        FissionYield {
            mass_number: 135,
            yield_fraction: 0.0650
        },
        FissionYield {
            mass_number: 137,
            yield_fraction: 0.0630
        },
        FissionYield {
            mass_number: 140,
            yield_fraction: 0.0620
        },
        FissionYield {
            mass_number: 144,
            yield_fraction: 0.0540
        },
        FissionYield {
            mass_number: 147,
            yield_fraction: 0.0224
        },
        FissionYield {
            mass_number: 151,
            yield_fraction: 0.0042
        },
        FissionYield {
            mass_number: 155,
            yield_fraction: 0.0003
        },
    ]
}

/// Returns thermal neutron fission yield distribution for Pu-239.
///
/// Source: ENDF/B-VIII.0 fission yield data.
#[must_use]
pub fn pu239_fission_yields() -> alloc::vec::Vec<FissionYield> {
    alloc::vec![
        FissionYield {
            mass_number: 85,
            yield_fraction: 0.0054
        },
        FissionYield {
            mass_number: 90,
            yield_fraction: 0.0211
        },
        FissionYield {
            mass_number: 95,
            yield_fraction: 0.0484
        },
        FissionYield {
            mass_number: 99,
            yield_fraction: 0.0620
        },
        FissionYield {
            mass_number: 103,
            yield_fraction: 0.0700
        },
        FissionYield {
            mass_number: 106,
            yield_fraction: 0.0425
        },
        FissionYield {
            mass_number: 110,
            yield_fraction: 0.0040
        },
        FissionYield {
            mass_number: 131,
            yield_fraction: 0.0385
        },
        FissionYield {
            mass_number: 134,
            yield_fraction: 0.0708
        },
        FissionYield {
            mass_number: 137,
            yield_fraction: 0.0666
        },
        FissionYield {
            mass_number: 140,
            yield_fraction: 0.0538
        },
        FissionYield {
            mass_number: 144,
            yield_fraction: 0.0375
        },
        FissionYield {
            mass_number: 147,
            yield_fraction: 0.0210
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dt_fusion_q_value() {
        let rxn = dt_fusion();
        assert!(
            (rxn.q_value_mev - 17.6).abs() < 0.1,
            "DT Q={} MeV",
            rxn.q_value_mev
        );
    }

    #[test]
    fn dd_fusion_q_value() {
        let rxn = dd_fusion_he3();
        assert!(
            (rxn.q_value_mev - 3.27).abs() < 0.1,
            "DD Q={} MeV",
            rxn.q_value_mev
        );
    }

    #[test]
    fn u235_fission_q_value() {
        let rxn = u235_fission();
        assert!(
            (rxn.q_value_mev - 200.0).abs() < 10.0,
            "U-235 fission Q={} MeV",
            rxn.q_value_mev
        );
    }

    #[test]
    fn q_value_computation() {
        // Simple test: if reactants weigh more than products, Q > 0
        let q = q_value(&[100.0, 200.0], &[290.0]);
        assert!((q - 10.0).abs() < 1e-10);
    }

    #[test]
    fn coulomb_barrier_dt() {
        let d = Nucleus::new(1, 2).unwrap();
        let t = Nucleus::new(1, 3).unwrap();
        let v = coulomb_barrier(&d, &t).unwrap();
        // Coulomb barrier for D-T should be a few hundred keV
        assert!(v > 0.1, "D-T barrier={v} MeV too low");
        assert!(v < 2.0, "D-T barrier={v} MeV too high");
    }

    #[test]
    fn coulomb_barrier_increases_with_z() {
        let d = Nucleus::new(1, 2).unwrap();
        let he4 = Nucleus::helium_4();
        let c12 = Nucleus::carbon_12();

        let v_d_he = coulomb_barrier(&d, &he4).unwrap();
        let v_d_c = coulomb_barrier(&d, &c12).unwrap();
        assert!(v_d_c > v_d_he, "C barrier should be higher than He barrier");
    }

    #[test]
    fn pp_chain_q_value() {
        let rxn = pp_chain_step1();
        assert!(
            (rxn.q_value_mev - 0.42).abs() < 0.05,
            "pp chain Q={} MeV",
            rxn.q_value_mev
        );
    }

    #[test]
    fn triple_alpha_q_value() {
        let rxn = triple_alpha();
        assert!(
            (rxn.q_value_mev - 7.275).abs() < 0.1,
            "Triple-alpha Q={} MeV",
            rxn.q_value_mev
        );
    }

    #[test]
    fn serde_roundtrip_reaction() {
        let rxn = dt_fusion();
        let json = serde_json::to_string(&rxn).unwrap();
        let back: NuclearReaction = serde_json::from_str(&json).unwrap();
        assert_eq!(rxn.name, back.name);
        assert!((rxn.q_value_mev - back.q_value_mev).abs() < 1e-10);
    }

    #[test]
    fn cno_cycle_q_value() {
        let rxn = cno_cycle();
        assert!(
            (rxn.q_value_mev - 25.03).abs() < 0.5,
            "CNO Q={} MeV",
            rxn.q_value_mev
        );
    }

    #[test]
    fn coulomb_constant_consistency() {
        use crate::constants::{COULOMB_MEV_FM, FINE_STRUCTURE};
        let hbar_c = 197.326_980_4;
        let expected = FINE_STRUCTURE * hbar_c;
        let rel_err = (COULOMB_MEV_FM - expected).abs() / expected;
        assert!(rel_err < 1e-5);
    }

    // --- Cross-section tests ---

    #[test]
    fn geometric_cross_section_fe56() {
        let fe56 = Nucleus::iron_56();
        let sigma = geometric_cross_section_barns(&fe56);
        // R ≈ 4.59 fm, σ = π R² ≈ 66.2 fm² ≈ 0.662 barns
        assert!(sigma > 0.5, "Fe-56 σ_geo={sigma} b too low");
        assert!(sigma < 1.0, "Fe-56 σ_geo={sigma} b too high");
    }

    #[test]
    fn geometric_cross_section_scales_with_a() {
        let c12 = Nucleus::carbon_12();
        let u238 = Nucleus::uranium_238();
        assert!(
            geometric_cross_section_barns(&u238) > geometric_cross_section_barns(&c12),
            "U-238 should have larger cross-section than C-12"
        );
    }

    #[test]
    fn thermal_u235_fission_cross_section() {
        let xs = thermal_neutron_cross_sections();
        let u235 = xs
            .iter()
            .find(|x| x.target == Nucleus::uranium_235())
            .unwrap();
        assert!(
            (u235.fission_barns - 585.0).abs() < 5.0,
            "U-235 σ_f={} b",
            u235.fission_barns
        );
    }

    #[test]
    fn serde_roundtrip_thermal_cross_section() {
        let xs = &thermal_neutron_cross_sections()[0];
        let json = serde_json::to_string(xs).unwrap();
        let back: ThermalCrossSection = serde_json::from_str(&json).unwrap();
        assert_eq!(xs.target, back.target);
        assert!((xs.fission_barns - back.fission_barns).abs() < 1e-10);
    }

    // --- Nucleosynthesis tests ---

    #[test]
    fn s_process_ends_at_pb208() {
        let pathway = s_process_main();
        assert_eq!(pathway.process, NucleosynthesisProcess::SProcess);
        let last = pathway.steps.last().unwrap();
        assert_eq!(last.nucleus.z(), 82, "s-process should end at Pb");
        assert_eq!(last.nucleus.a(), 208, "s-process should end at Pb-208");
    }

    #[test]
    fn r_process_produces_actinides() {
        let pathway = r_process_main();
        assert_eq!(pathway.process, NucleosynthesisProcess::RProcess);
        let last = pathway.steps.last().unwrap();
        assert_eq!(last.nucleus.z(), 92, "r-process should reach U");
    }

    #[test]
    fn s_process_peaks_at_magic_neutrons() {
        let pathway = s_process_main();
        // Sr-88 has N=50 (magic), Ba-138 has N=82 (magic), Pb-208 has N=126 (magic)
        let sr = pathway.steps.iter().find(|s| s.nucleus.a() == 88).unwrap();
        assert_eq!(sr.nucleus.n(), 50, "Sr-88 should have N=50");
        let ba = pathway.steps.iter().find(|s| s.nucleus.a() == 138).unwrap();
        assert_eq!(ba.nucleus.n(), 82, "Ba-138 should have N=82");
    }

    #[test]
    fn serde_roundtrip_nucleosynthesis() {
        let pathway = s_process_main();
        let json = serde_json::to_string(&pathway).unwrap();
        let back: NucleosynthesisPathway = serde_json::from_str(&json).unwrap();
        assert_eq!(pathway.name, back.name);
        assert_eq!(pathway.steps.len(), back.steps.len());
    }

    // --- Neutron moderation tests ---

    #[test]
    fn hydrogen_max_energy_loss() {
        // Hydrogen (A=1): max energy loss = 100% (head-on billiard)
        let frac = max_energy_loss_fraction(1);
        assert!((frac - 1.0).abs() < 1e-10);
    }

    #[test]
    fn heavy_nucleus_small_energy_loss() {
        // U-238: max loss ≈ 4/A ≈ 1.7%
        let frac = max_energy_loss_fraction(238);
        assert!(frac < 0.02, "U-238 max loss={frac}");
        assert!(frac > 0.01, "U-238 max loss={frac}");
    }

    #[test]
    fn hydrogen_lethargy_gain() {
        let xi = average_lethargy_gain(1);
        assert!((xi - 1.0).abs() < 1e-10, "H ξ={xi}");
    }

    #[test]
    fn collisions_hydrogen_thermalize() {
        // 2 MeV -> 0.0253 eV in hydrogen: ~18 collisions
        let nc = collisions_to_thermalize(1, 2e6, 0.0253);
        assert!(nc > 15.0 && nc < 25.0, "H collisions={nc}");
    }

    #[test]
    fn collisions_graphite_more_than_hydrogen() {
        let nc_h = collisions_to_thermalize(1, 2e6, 0.0253);
        let nc_c = collisions_to_thermalize(12, 2e6, 0.0253);
        assert!(nc_c > nc_h, "C-12 needs more collisions than H");
    }

    // --- Fission yield tests ---

    #[test]
    fn u235_fission_yield_double_hump() {
        let yields = u235_fission_yields();
        // Light peak around A=95, heavy peak around A=137
        let light_peak = yields.iter().find(|y| y.mass_number == 95).unwrap();
        let valley = yields.iter().find(|y| y.mass_number == 115).unwrap();
        let heavy_peak = yields.iter().find(|y| y.mass_number == 137).unwrap();
        assert!(light_peak.yield_fraction > valley.yield_fraction);
        assert!(heavy_peak.yield_fraction > valley.yield_fraction);
    }

    #[test]
    fn pu239_fission_yields_exist() {
        let yields = pu239_fission_yields();
        assert!(yields.len() >= 10);
    }

    #[test]
    fn serde_roundtrip_fission_yield() {
        let fy = &u235_fission_yields()[0];
        let json = serde_json::to_string(fy).unwrap();
        let back: FissionYield = serde_json::from_str(&json).unwrap();
        assert_eq!(fy.mass_number, back.mass_number);
    }
}
