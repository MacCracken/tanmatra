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

    let y = seconds_per_year;
    let d = seconds_per_day;
    let m = seconds_per_minute;

    alloc::vec![
        gs(1, 3, "H-3", 12.32 * y, DecayMode::BetaMinus),
        gs(6, 14, "C-14", 5730.0 * y, DecayMode::BetaMinus),
        gs(19, 40, "K-40", 1.248e9 * y, DecayMode::BetaMinus),
        gs(27, 60, "Co-60", 5.2714 * y, DecayMode::BetaMinus),
        gs(38, 90, "Sr-90", 28.79 * y, DecayMode::BetaMinus),
        gs(53, 131, "I-131", 8.0252 * d, DecayMode::BetaMinus),
        gs(55, 137, "Cs-137", 30.05 * y, DecayMode::BetaMinus),
        gs(86, 222, "Rn-222", 3.8222 * d, DecayMode::Alpha),
        gs(88, 226, "Ra-226", 1600.0 * y, DecayMode::Alpha),
        gs(90, 230, "Th-230", 75_380.0 * y, DecayMode::Alpha),
        gs(90, 234, "Th-234", 24.10 * d, DecayMode::BetaMinus),
        // Pa-234m: isomeric state, excitation energy 73.92 keV (NNDC)
        iso(91, 234, "Pa-234m", 1.17 * m, DecayMode::BetaMinus, 73.92),
        gs(92, 234, "U-234", 245_250.0 * y, DecayMode::Alpha),
        gs(92, 235, "U-235", 7.04e8 * y, DecayMode::Alpha),
        gs(92, 238, "U-238", 4.468e9 * y, DecayMode::Alpha),
        // --- Full U-238 decay chain intermediates (NNDC) ---
        gs(84, 218, "Po-218", 3.098 * m, DecayMode::Alpha),
        gs(82, 214, "Pb-214", 26.8 * m, DecayMode::BetaMinus),
        gs(83, 214, "Bi-214", 19.9 * m, DecayMode::BetaMinus),
        gs(84, 214, "Po-214", 164.3e-6, DecayMode::Alpha), // 164.3 µs
        gs(82, 210, "Pb-210", 22.2 * y, DecayMode::BetaMinus),
        gs(83, 210, "Bi-210", 5.012 * d, DecayMode::BetaMinus),
        gs(84, 210, "Po-210", 138.376 * d, DecayMode::Alpha),
        // Pb-206 is stable (end of U-238 chain) — not included
        // --- Other notable isotopes ---
        gs(94, 239, "Pu-239", 24_110.0 * y, DecayMode::Alpha),
        gs(95, 241, "Am-241", 432.6 * y, DecayMode::Alpha),
        // Notable isomers
        // Ta-180m: longest-lived nuclear isomer (>1.2e15 years, effectively stable)
        iso(
            73,
            180,
            "Ta-180m",
            1.2e15 * y,
            DecayMode::IsomericTransition,
            77.1
        ),
        // Hf-178m2: high-spin isomer
        iso(
            72,
            178,
            "Hf-178m2",
            31.0 * y,
            DecayMode::IsomericTransition,
            2446.1
        ),
        // Tc-99m: medical imaging isotope
        iso(
            43,
            99,
            "Tc-99m",
            6.006 * 3600.0,
            DecayMode::IsomericTransition,
            142.68
        ),
        // Am-242m: used in nuclear reactors
        iso(
            95,
            242,
            "Am-242m",
            141.0 * y,
            DecayMode::IsomericTransition,
            48.6
        ),
    ]
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
    fn serde_roundtrip_isomeric_transition() {
        let mode = DecayMode::IsomericTransition;
        let json = serde_json::to_string(&mode).unwrap();
        let back: DecayMode = serde_json::from_str(&json).unwrap();
        assert_eq!(mode, back);
    }
}
