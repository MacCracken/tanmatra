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
        // COULOMB_MEV_FM = alpha * hbar * c ≈ (1/137.036) * 197.327 ≈ 1.4398 MeV*fm
        use crate::constants::{COULOMB_MEV_FM, FINE_STRUCTURE};
        let hbar_c = 197.326_980_4; // MeV*fm (CODATA 2022)
        let expected = FINE_STRUCTURE * hbar_c;
        let rel_err = (COULOMB_MEV_FM - expected).abs() / expected;
        assert!(
            rel_err < 1e-5,
            "COULOMB_MEV_FM={COULOMB_MEV_FM} vs derived={expected}"
        );
    }
}
