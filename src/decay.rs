//! Radioactive decay calculations.
//!
//! Provides decay modes, decay constant calculations, activity computations,
//! and a library of known isotopes with real half-lives from NNDC/NUBASE.

use crate::error::TanmatraError;
use crate::nucleus::Nucleus;
use alloc::string::String;
use alloc::vec::Vec;
use serde::{Deserialize, Serialize};

/// Natural logarithm of 2.
const LN2: f64 = core::f64::consts::LN_2;

/// Radioactive decay modes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum DecayMode {
    /// Alpha decay: emits He-4 nucleus (Z-2, A-4).
    Alpha,
    /// Beta-minus decay: neutron -> proton + electron + antineutrino (Z+1, A).
    BetaMinus,
    /// Beta-plus decay: proton -> neutron + positron + neutrino (Z-1, A).
    BetaPlus,
    /// Electron capture: proton + electron -> neutron + neutrino (Z-1, A).
    ElectronCapture,
    /// Gamma decay: excited nucleus emits photon (Z, A unchanged).
    Gamma,
    /// Spontaneous fission.
    SpontaneousFission,
    /// Proton emission (Z-1, A-1).
    ProtonEmission,
    /// Neutron emission (Z, A-1).
    NeutronEmission,
    /// Isomeric transition: excited state -> lower state via gamma emission.
    IsomericTransition,
}

/// A known radioactive isotope with its half-life and primary decay mode.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Isotope {
    /// The nucleus.
    pub nucleus: Nucleus,
    /// Name/symbol of the isotope.
    pub name: String,
    /// Half-life in seconds.
    pub half_life_seconds: f64,
    /// Primary decay mode.
    pub primary_decay: DecayMode,
    /// Whether this is an isomeric (metastable) state.
    pub is_isomer: bool,
    /// Excitation energy above ground state in keV (0.0 for ground state).
    pub excitation_energy_kev: f64,
}

/// Returns the decay constant lambda = ln(2) / t_half.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidHalfLife`] if `half_life` is not positive and finite.
#[inline]
pub fn decay_constant(half_life: f64) -> Result<f64, TanmatraError> {
    if half_life <= 0.0 || !half_life.is_finite() {
        return Err(TanmatraError::InvalidHalfLife(alloc::format!(
            "{half_life} seconds is not a valid half-life"
        )));
    }
    Ok(LN2 / half_life)
}

/// Returns the fraction of material remaining after time `t` given `half_life`.
///
/// N(t)/N0 = (1/2)^(t/t_half) = exp(-lambda * t)
///
/// Returns 0.0 if `half_life` is not positive/finite or `time` is negative.
#[must_use]
#[inline]
pub fn remaining_fraction(half_life: f64, time: f64) -> f64 {
    if half_life <= 0.0 || !half_life.is_finite() || time < 0.0 {
        return 0.0;
    }
    libm::pow(0.5, time / half_life)
}

/// Returns the activity in becquerels (decays per second).
///
/// A = lambda * N = (ln2 / t_half) * N
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidHalfLife`] if `half_life` is not positive and finite.
#[inline]
pub fn activity_bq(half_life: f64, num_atoms: f64) -> Result<f64, TanmatraError> {
    let lambda = decay_constant(half_life)?;
    Ok(lambda * num_atoms)
}

/// Performs alpha decay on a nucleus, returning the daughter nucleus.
///
/// Alpha decay: (Z, A) -> (Z-2, A-4) + He-4
///
/// # Errors
///
/// Returns [`TanmatraError::DecayNotPossible`] if Z < 3 or A < 5.
pub fn alpha_decay(parent: &Nucleus) -> Result<Nucleus, TanmatraError> {
    if parent.z() < 3 || parent.a() < 5 {
        return Err(TanmatraError::DecayNotPossible(alloc::format!(
            "alpha decay requires Z >= 3 and A >= 5, got Z={} A={}",
            parent.z(),
            parent.a()
        )));
    }
    // Safe: validated above
    Nucleus::new(parent.z() - 2, parent.a() - 4)
        .map_err(|e| TanmatraError::DecayNotPossible(alloc::format!("{e}")))
}

/// Performs beta-minus decay on a nucleus, returning the daughter nucleus.
///
/// Beta-minus: (Z, A) -> (Z+1, A) + e- + antineutrino
///
/// # Errors
///
/// Returns [`TanmatraError::DecayNotPossible`] if the nucleus has no neutrons.
pub fn beta_minus_decay(parent: &Nucleus) -> Result<Nucleus, TanmatraError> {
    if parent.n() == 0 {
        return Err(TanmatraError::DecayNotPossible(
            alloc::string::String::from("beta-minus decay requires neutrons, nucleus has N=0"),
        ));
    }
    Nucleus::new(parent.z() + 1, parent.a())
        .map_err(|e| TanmatraError::DecayNotPossible(alloc::format!("{e}")))
}

/// Performs beta-plus decay on a nucleus, returning the daughter nucleus.
///
/// Beta-plus: (Z, A) -> (Z-1, A) + e+ + neutrino
///
/// # Errors
///
/// Returns [`TanmatraError::DecayNotPossible`] if Z < 2.
pub fn beta_plus_decay(parent: &Nucleus) -> Result<Nucleus, TanmatraError> {
    if parent.z() < 2 {
        return Err(TanmatraError::DecayNotPossible(alloc::format!(
            "beta-plus decay requires Z >= 2, got Z={}",
            parent.z()
        )));
    }
    Nucleus::new(parent.z() - 1, parent.a())
        .map_err(|e| TanmatraError::DecayNotPossible(alloc::format!("{e}")))
}

/// Generates a decay chain from the given isotope until a stable nucleus
/// is reached or `max_steps` is exceeded.
///
/// Uses the known isotope database to determine decay modes and daughters.
/// Returns a list of (Nucleus, DecayMode) pairs for each step.
#[must_use]
pub fn decay_chain(start: &Nucleus, max_steps: usize) -> Vec<(Nucleus, DecayMode)> {
    let mut chain = Vec::new();
    let mut current = *start;
    let all = known_isotopes();

    for _ in 0..max_steps {
        // Look up the isotope in our known database
        let isotope = all.iter().find(|iso| iso.nucleus == current);

        match isotope {
            Some(iso) => {
                let mode = iso.primary_decay;
                let daughter = match mode {
                    DecayMode::Alpha => alpha_decay(&current),
                    DecayMode::BetaMinus => beta_minus_decay(&current),
                    DecayMode::BetaPlus | DecayMode::ElectronCapture => beta_plus_decay(&current),
                    DecayMode::IsomericTransition | DecayMode::Gamma => {
                        // IT/gamma: nucleus unchanged, transition to ground state
                        Ok(current)
                    }
                    _ => break, // Fission, proton/neutron emission -- stop chain
                };

                match daughter {
                    Ok(d) => {
                        chain.push((current, mode));
                        current = d;
                    }
                    Err(_) => break,
                }
            }
            None => break, // Unknown or stable
        }
    }

    chain
}

/// Returns a collection of known radioactive isotopes with real half-lives.
///
/// Half-lives are from the NNDC (National Nuclear Data Center) / NUBASE evaluations.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn known_isotopes() -> Vec<Isotope> {
    let seconds_per_year = 365.25 * 24.0 * 3600.0;
    let seconds_per_day = 24.0 * 3600.0;
    let seconds_per_minute = 60.0;

    // Helper to create a ground-state isotope
    let gs = |z: u32, a: u32, name: &str, half_life: f64, mode: DecayMode| Isotope {
        nucleus: Nucleus::new(z, a).unwrap_or_else(|_| Nucleus::hydrogen_1()),
        name: String::from(name),
        half_life_seconds: half_life,
        primary_decay: mode,
        is_isomer: false,
        excitation_energy_kev: 0.0,
    };
    let iso = |z: u32, a: u32, name: &str, half_life: f64, mode: DecayMode, exc_kev: f64| Isotope {
        nucleus: Nucleus::new(z, a).unwrap_or_else(|_| Nucleus::hydrogen_1()),
        name: String::from(name),
        half_life_seconds: half_life,
        primary_decay: mode,
        is_isomer: true,
        excitation_energy_kev: exc_kev,
    };

    let yr = seconds_per_year;
    let dy = seconds_per_day;
    let mn = seconds_per_minute;

    let hr = 3600.0; // seconds per hour
    let sc = 1.0; // seconds

    alloc::vec![
        // === Cosmogenic / Light ===
        gs(1, 3, "H-3", 12.32 * yr, DecayMode::BetaMinus),
        gs(4, 7, "Be-7", 53.22 * dy, DecayMode::ElectronCapture),
        gs(4, 10, "Be-10", 1.387e6 * yr, DecayMode::BetaMinus),
        gs(6, 14, "C-14", 5730.0 * yr, DecayMode::BetaMinus),
        gs(9, 18, "F-18", 109.77 * mn, DecayMode::BetaPlus),
        gs(11, 22, "Na-22", 2.6029 * yr, DecayMode::BetaPlus),
        gs(13, 26, "Al-26", 7.17e5 * yr, DecayMode::BetaPlus),
        gs(14, 32, "Si-32", 153.0 * yr, DecayMode::BetaMinus),
        gs(15, 32, "P-32", 14.268 * dy, DecayMode::BetaMinus),
        gs(16, 35, "S-35", 87.37 * dy, DecayMode::BetaMinus),
        gs(17, 36, "Cl-36", 3.01e5 * yr, DecayMode::BetaMinus),
        gs(18, 39, "Ar-39", 268.0 * yr, DecayMode::BetaMinus),
        gs(19, 40, "K-40", 1.248e9 * yr, DecayMode::BetaMinus),
        gs(20, 41, "Ca-41", 9.94e4 * yr, DecayMode::ElectronCapture),
        gs(24, 51, "Cr-51", 27.701 * dy, DecayMode::ElectronCapture),
        gs(25, 53, "Mn-53", 3.74e6 * yr, DecayMode::ElectronCapture),
        gs(26, 55, "Fe-55", 2.744 * yr, DecayMode::ElectronCapture),
        gs(26, 59, "Fe-59", 44.490 * dy, DecayMode::BetaMinus),
        gs(27, 57, "Co-57", 271.74 * dy, DecayMode::ElectronCapture),
        gs(27, 60, "Co-60", 5.2714 * yr, DecayMode::BetaMinus),
        gs(28, 63, "Ni-63", 101.2 * yr, DecayMode::BetaMinus),
        gs(29, 64, "Cu-64", 12.701 * hr, DecayMode::BetaPlus),
        gs(30, 65, "Zn-65", 243.93 * dy, DecayMode::ElectronCapture),
        gs(31, 67, "Ga-67", 3.2617 * dy, DecayMode::ElectronCapture),
        gs(31, 68, "Ga-68", 67.71 * mn, DecayMode::BetaPlus),
        gs(36, 85, "Kr-85", 10.739 * yr, DecayMode::BetaMinus),
        gs(37, 87, "Rb-87", 4.923e10 * yr, DecayMode::BetaMinus),
        gs(38, 89, "Sr-89", 50.563 * dy, DecayMode::BetaMinus),
        gs(38, 90, "Sr-90", 28.79 * yr, DecayMode::BetaMinus),
        gs(39, 90, "Y-90", 64.053 * hr, DecayMode::BetaMinus),
        gs(40, 95, "Zr-95", 64.032 * dy, DecayMode::BetaMinus),
        gs(41, 95, "Nb-95", 34.991 * dy, DecayMode::BetaMinus),
        gs(42, 99, "Mo-99", 65.924 * hr, DecayMode::BetaMinus),
        gs(43, 99, "Tc-99", 2.111e5 * yr, DecayMode::BetaMinus),
        gs(44, 103, "Ru-103", 39.247 * dy, DecayMode::BetaMinus),
        gs(44, 106, "Ru-106", 371.8 * dy, DecayMode::BetaMinus),
        gs(48, 109, "Cd-109", 461.9 * dy, DecayMode::ElectronCapture),
        gs(49, 111, "In-111", 2.8047 * dy, DecayMode::ElectronCapture),
        gs(50, 113, "Sn-113", 115.09 * dy, DecayMode::ElectronCapture),
        gs(51, 125, "Sb-125", 2.7586 * yr, DecayMode::BetaMinus),
        gs(52, 132, "Te-132", 3.204 * dy, DecayMode::BetaMinus),
        gs(53, 123, "I-123", 13.2235 * hr, DecayMode::ElectronCapture),
        gs(53, 125, "I-125", 59.407 * dy, DecayMode::ElectronCapture),
        gs(53, 129, "I-129", 1.57e7 * yr, DecayMode::BetaMinus),
        gs(53, 131, "I-131", 8.0252 * dy, DecayMode::BetaMinus),
        gs(55, 134, "Cs-134", 2.0652 * yr, DecayMode::BetaMinus),
        gs(55, 137, "Cs-137", 30.05 * yr, DecayMode::BetaMinus),
        gs(56, 133, "Ba-133", 10.551 * yr, DecayMode::ElectronCapture),
        gs(56, 140, "Ba-140", 12.752 * dy, DecayMode::BetaMinus),
        gs(57, 140, "La-140", 1.6781 * dy, DecayMode::BetaMinus),
        gs(58, 144, "Ce-144", 284.91 * dy, DecayMode::BetaMinus),
        gs(61, 147, "Pm-147", 2.6234 * yr, DecayMode::BetaMinus),
        gs(62, 153, "Sm-153", 46.284 * hr, DecayMode::BetaMinus),
        gs(63, 152, "Eu-152", 13.517 * yr, DecayMode::ElectronCapture),
        gs(63, 154, "Eu-154", 8.601 * yr, DecayMode::BetaMinus),
        gs(63, 155, "Eu-155", 4.7611 * yr, DecayMode::BetaMinus),
        gs(71, 177, "Lu-177", 6.6463 * dy, DecayMode::BetaMinus),
        gs(75, 186, "Re-186", 3.7186 * dy, DecayMode::BetaMinus),
        gs(75, 188, "Re-188", 17.004 * hr, DecayMode::BetaMinus),
        gs(77, 192, "Ir-192", 73.829 * dy, DecayMode::BetaMinus),
        gs(79, 198, "Au-198", 2.6941 * dy, DecayMode::BetaMinus),
        gs(81, 204, "Tl-204", 3.783 * yr, DecayMode::BetaMinus),
        // === Th-232 decay chain ===
        gs(90, 232, "Th-232", 1.405e10 * yr, DecayMode::Alpha),
        gs(88, 228, "Ra-228", 5.75 * yr, DecayMode::BetaMinus),
        gs(89, 228, "Ac-228", 6.15 * hr, DecayMode::BetaMinus),
        gs(90, 228, "Th-228", 1.9125 * yr, DecayMode::Alpha),
        gs(88, 224, "Ra-224", 3.6319 * dy, DecayMode::Alpha),
        gs(86, 220, "Rn-220", 55.6 * sc, DecayMode::Alpha),
        gs(84, 216, "Po-216", 0.145 * sc, DecayMode::Alpha),
        gs(82, 212, "Pb-212", 10.64 * hr, DecayMode::BetaMinus),
        gs(83, 212, "Bi-212", 60.55 * mn, DecayMode::Alpha),
        gs(81, 208, "Tl-208", 3.053 * mn, DecayMode::BetaMinus),
        gs(84, 212, "Po-212", 0.299e-6, DecayMode::Alpha),
        // === U-235 (actinium) decay chain ===
        gs(90, 231, "Th-231", 25.52 * hr, DecayMode::BetaMinus),
        gs(91, 231, "Pa-231", 3.276e4 * yr, DecayMode::Alpha),
        gs(89, 227, "Ac-227", 21.772 * yr, DecayMode::BetaMinus),
        gs(90, 227, "Th-227", 18.697 * dy, DecayMode::Alpha),
        gs(87, 223, "Fr-223", 22.00 * mn, DecayMode::BetaMinus),
        gs(88, 223, "Ra-223", 11.43 * dy, DecayMode::Alpha),
        gs(86, 219, "Rn-219", 3.96 * sc, DecayMode::Alpha),
        gs(84, 215, "Po-215", 1.781e-3, DecayMode::Alpha),
        gs(82, 211, "Pb-211", 36.1 * mn, DecayMode::BetaMinus),
        gs(83, 211, "Bi-211", 2.14 * mn, DecayMode::Alpha),
        gs(81, 207, "Tl-207", 4.77 * mn, DecayMode::BetaMinus),
        // === U-238 decay chain ===
        gs(92, 238, "U-238", 4.468e9 * yr, DecayMode::Alpha),
        gs(90, 234, "Th-234", 24.10 * dy, DecayMode::BetaMinus),
        iso(91, 234, "Pa-234m", 1.17 * mn, DecayMode::BetaMinus, 73.92),
        gs(91, 234, "Pa-234", 6.70 * hr, DecayMode::BetaMinus),
        gs(92, 234, "U-234", 245_250.0 * yr, DecayMode::Alpha),
        gs(92, 235, "U-235", 7.04e8 * yr, DecayMode::Alpha),
        gs(90, 230, "Th-230", 75_380.0 * yr, DecayMode::Alpha),
        gs(88, 226, "Ra-226", 1600.0 * yr, DecayMode::Alpha),
        gs(86, 222, "Rn-222", 3.8222 * dy, DecayMode::Alpha),
        gs(84, 218, "Po-218", 3.098 * mn, DecayMode::Alpha),
        gs(82, 214, "Pb-214", 26.8 * mn, DecayMode::BetaMinus),
        gs(83, 214, "Bi-214", 19.9 * mn, DecayMode::BetaMinus),
        gs(84, 214, "Po-214", 164.3e-6, DecayMode::Alpha),
        gs(82, 210, "Pb-210", 22.2 * yr, DecayMode::BetaMinus),
        gs(83, 210, "Bi-210", 5.012 * dy, DecayMode::BetaMinus),
        gs(84, 210, "Po-210", 138.376 * dy, DecayMode::Alpha),
        // === Actinides ===
        gs(89, 225, "Ac-225", 9.920 * dy, DecayMode::Alpha),
        gs(93, 237, "Np-237", 2.144e6 * yr, DecayMode::Alpha),
        gs(94, 238, "Pu-238", 87.74 * yr, DecayMode::Alpha),
        gs(94, 239, "Pu-239", 24_110.0 * yr, DecayMode::Alpha),
        gs(94, 240, "Pu-240", 6_561.0 * yr, DecayMode::Alpha),
        gs(94, 241, "Pu-241", 14.329 * yr, DecayMode::BetaMinus),
        gs(94, 242, "Pu-242", 3.75e5 * yr, DecayMode::Alpha),
        gs(95, 241, "Am-241", 432.6 * yr, DecayMode::Alpha),
        gs(96, 244, "Cm-244", 18.11 * yr, DecayMode::Alpha),
        gs(98, 252, "Cf-252", 2.645 * yr, DecayMode::Alpha),
        // === Isomers ===
        iso(
            43,
            99,
            "Tc-99m",
            6.006 * hr,
            DecayMode::IsomericTransition,
            142.68
        ),
        iso(
            72,
            178,
            "Hf-178m2",
            31.0 * yr,
            DecayMode::IsomericTransition,
            2446.1
        ),
        iso(
            73,
            180,
            "Ta-180m",
            1.2e15 * yr,
            DecayMode::IsomericTransition,
            77.1
        ),
        iso(
            95,
            242,
            "Am-242m",
            141.0 * yr,
            DecayMode::IsomericTransition,
            48.6
        ),
    ]
}

// ---------------------------------------------------------------------------
// Bateman equations for sequential decay chains
// ---------------------------------------------------------------------------

/// Computes nuclide populations in a sequential decay chain at time t
/// using the Bateman equations.
///
/// For a chain: N₁ → N₂ → N₃ → ... → N_n (stable)
///
/// Given initial atoms N₁(0) and decay constants λ₁, λ₂, ..., λ_{n-1},
/// the Bateman solution gives the number of atoms of each species at time t.
///
/// Parameters:
/// - `decay_constants`: λ values in s⁻¹ for each species (last entry = 0 means stable)
/// - `initial_atoms`: N₁(0) (only the first species has nonzero initial population)
/// - `time_seconds`: time at which to evaluate
///
/// Returns a Vec of atom counts for each species in the chain.
#[must_use]
#[allow(clippy::needless_range_loop)]
pub fn bateman_chain(decay_constants: &[f64], initial_atoms: f64, time_seconds: f64) -> Vec<f64> {
    let chain_len = decay_constants.len();
    if chain_len == 0 {
        return Vec::new();
    }

    let mut populations = alloc::vec![0.0; chain_len];

    // First species: simple exponential decay
    let lambda1 = decay_constants[0];
    if lambda1 <= 0.0 {
        // First species is stable
        populations[0] = initial_atoms;
        return populations;
    }
    populations[0] = initial_atoms * libm::exp(-lambda1 * time_seconds);

    // Bateman solution for species i (1-indexed, but 0-indexed in array):
    // N_i(t) = N_1(0) * Π_{j=1}^{i-1} λ_j * Σ_{j=1}^{i} [exp(-λ_j t) / Π_{k≠j} (λ_k - λ_j)]
    for idx in 1..chain_len {
        let mut sum = 0.0;

        for jj in 0..=idx {
            let lambda_j = decay_constants[jj];

            // Compute product of (λ_k - λ_j) for k ≠ j, k in 0..=idx
            let mut product = 1.0;
            let mut degenerate = false;
            for kk in 0..=idx {
                if kk == jj {
                    continue;
                }
                let diff = decay_constants[kk] - lambda_j;
                if diff.abs() < 1e-30 {
                    degenerate = true;
                    break;
                }
                product *= diff;
            }

            if degenerate || product.abs() < 1e-100 {
                continue;
            }

            sum += libm::exp(-lambda_j * time_seconds) / product;
        }

        // Multiply by product of all λ_j for j < idx, and N_1(0)
        let mut lambda_product = 1.0;
        for jj in 0..idx {
            lambda_product *= decay_constants[jj];
        }

        populations[idx] = initial_atoms * lambda_product * sum;

        // Clamp negative values from numerical noise
        if populations[idx] < 0.0 {
            populations[idx] = 0.0;
        }
    }

    populations
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn c14_half_life() {
        let c14 = known_isotopes()
            .into_iter()
            .find(|i| i.name == "C-14")
            .unwrap();
        let years = c14.half_life_seconds / (365.25 * 24.0 * 3600.0);
        assert!((years - 5730.0).abs() < 1.0, "C-14 half-life={years} years");
    }

    #[test]
    fn u238_half_life() {
        let u238 = known_isotopes()
            .into_iter()
            .find(|i| i.name == "U-238")
            .unwrap();
        let years = u238.half_life_seconds / (365.25 * 24.0 * 3600.0);
        assert!(
            (years - 4.468e9).abs() < 1e7,
            "U-238 half-life={years} years"
        );
    }

    #[test]
    fn decay_constant_c14() {
        let t_half = 5730.0 * 365.25 * 24.0 * 3600.0;
        let lambda = decay_constant(t_half).unwrap();
        // lambda ≈ 3.83e-12 /s
        assert!(lambda > 3.8e-12);
        assert!(lambda < 3.9e-12);
    }

    #[test]
    fn remaining_fraction_one_half_life() {
        let frac = remaining_fraction(100.0, 100.0);
        assert!((frac - 0.5).abs() < 1e-10);
    }

    #[test]
    fn remaining_fraction_two_half_lives() {
        let frac = remaining_fraction(100.0, 200.0);
        assert!((frac - 0.25).abs() < 1e-10);
    }

    #[test]
    fn activity_bq_basic() {
        let t_half = 1.0; // 1 second half-life
        let n = 1e6;
        let a = activity_bq(t_half, n).unwrap();
        let expected = LN2 * 1e6;
        assert!((a - expected).abs() < 1.0);
    }

    #[test]
    fn u238_alpha_decay_to_th234() {
        let u238 = Nucleus::uranium_238();
        let daughter = alpha_decay(&u238).unwrap();
        assert_eq!(daughter.z(), 90); // Thorium
        assert_eq!(daughter.a(), 234);
    }

    #[test]
    fn beta_minus_increases_z() {
        let n = Nucleus::new(6, 14).unwrap(); // C-14
        let daughter = beta_minus_decay(&n).unwrap();
        assert_eq!(daughter.z(), 7); // N-14
        assert_eq!(daughter.a(), 14);
    }

    #[test]
    fn alpha_decay_h1_fails() {
        let h1 = Nucleus::hydrogen_1();
        assert!(alpha_decay(&h1).is_err());
    }

    #[test]
    fn invalid_half_life() {
        assert!(decay_constant(0.0).is_err());
        assert!(decay_constant(-1.0).is_err());
        assert!(decay_constant(f64::INFINITY).is_err());
        assert!(decay_constant(f64::NAN).is_err());
    }

    #[test]
    fn decay_chain_u238() {
        let u238 = Nucleus::uranium_238();
        let chain = decay_chain(&u238, 5);
        // U-238 -> Th-234 (alpha) -> Pa-234m (beta-) -> U-234 (beta-) -> Th-230 (alpha)
        assert!(!chain.is_empty());
        // First step should be alpha decay
        assert_eq!(chain[0].1, DecayMode::Alpha);
        // Daughter of first step should be Th-234
        if chain.len() > 1 {
            assert_eq!(chain[1].0.z(), 90);
            assert_eq!(chain[1].0.a(), 234);
        }
    }

    #[test]
    fn serde_roundtrip_decay_mode() {
        let mode = DecayMode::Alpha;
        let json = serde_json::to_string(&mode).unwrap();
        let back: DecayMode = serde_json::from_str(&json).unwrap();
        assert_eq!(mode, back);
    }

    #[test]
    fn serde_roundtrip_isotope() {
        let iso = &known_isotopes()[0];
        let json = serde_json::to_string(iso).unwrap();
        let back: Isotope = serde_json::from_str(&json).unwrap();
        assert_eq!(iso.name, back.name);
        assert!((iso.half_life_seconds - back.half_life_seconds).abs() < 1.0);
    }

    #[test]
    fn full_u238_chain_to_pb206() {
        let u238 = Nucleus::uranium_238();
        let chain = decay_chain(&u238, 20);
        // Full chain: U-238 -> Th-234 -> Pa-234m -> U-234 -> Th-230 -> Ra-226
        // -> Rn-222 -> Po-218 -> Pb-214 -> Bi-214 -> Po-214 -> Pb-210
        // -> Bi-210 -> Po-210 -> Pb-206 (stable, not in chain)
        // That's 14 steps
        assert!(
            chain.len() >= 14,
            "U-238 chain should have at least 14 steps, got {}",
            chain.len()
        );
        // Last nucleus in chain should decay to Pb-206
        if let Some(last) = chain.last() {
            let daughter = match last.1 {
                DecayMode::Alpha => alpha_decay(&last.0).ok(),
                DecayMode::BetaMinus => beta_minus_decay(&last.0).ok(),
                _ => None,
            };
            if let Some(d) = daughter {
                assert_eq!(d.z(), 82, "Chain should end at Pb (Z=82)");
                assert_eq!(d.a(), 206, "Chain should end at A=206");
            }
        }
    }

    #[test]
    fn tc99m_is_isomer() {
        let isotopes = known_isotopes();
        let tc99m = isotopes.iter().find(|i| i.name == "Tc-99m").unwrap();
        assert!(tc99m.is_isomer);
        assert!(tc99m.excitation_energy_kev > 140.0);
        assert_eq!(tc99m.primary_decay, DecayMode::IsomericTransition);
    }

    #[test]
    fn isotope_database_has_100_plus() {
        let isotopes = known_isotopes();
        assert!(
            isotopes.len() >= 100,
            "Expected 100+ isotopes, got {}",
            isotopes.len()
        );
    }

    #[test]
    fn serde_roundtrip_isomeric_transition() {
        let mode = DecayMode::IsomericTransition;
        let json = serde_json::to_string(&mode).unwrap();
        let back: DecayMode = serde_json::from_str(&json).unwrap();
        assert_eq!(mode, back);
    }

    // --- Bateman equation tests ---

    #[test]
    fn bateman_single_species_decay() {
        // Single species with λ=ln2 (t_half = 1s), after 1 half-life
        let lambdas = [LN2];
        let pops = bateman_chain(&lambdas, 1000.0, 1.0);
        assert_eq!(pops.len(), 1);
        assert!((pops[0] - 500.0).abs() < 1.0, "N={}, expected 500", pops[0]);
    }

    #[test]
    fn bateman_two_species_chain() {
        // A -> B (stable): λ_A = ln2, λ_B = 0
        // After t=1s: N_A = 500, N_B = 500
        let lambdas = [LN2, 0.0];
        let pops = bateman_chain(&lambdas, 1000.0, 1.0);
        assert_eq!(pops.len(), 2);
        assert!((pops[0] - 500.0).abs() < 1.0);
        assert!((pops[1] - 500.0).abs() < 5.0, "N_B={}", pops[1]);
    }

    #[test]
    fn bateman_conservation() {
        // Total atoms should be conserved (last species is stable)
        let lambdas = [0.1, 0.05, 0.0]; // A -> B -> C (stable)
        let n0 = 1e6;
        let pops = bateman_chain(&lambdas, n0, 100.0);
        let total: f64 = pops.iter().sum();
        let rel_err = (total - n0).abs() / n0;
        assert!(rel_err < 0.01, "Total={total}, expected {n0}");
    }

    #[test]
    fn bateman_empty_chain() {
        let pops = bateman_chain(&[], 1000.0, 1.0);
        assert!(pops.is_empty());
    }
}
