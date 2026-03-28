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
                    _ => break, // Gamma, fission, etc. -- stop chain
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

    alloc::vec![
        // Tritium: H-3, t½ = 12.32 years
        Isotope {
            nucleus: Nucleus::new(1, 3).unwrap_or_else(|_| Nucleus::hydrogen_1()),
            name: String::from("H-3"),
            half_life_seconds: 12.32 * seconds_per_year,
            primary_decay: DecayMode::BetaMinus,
        },
        // Carbon-14: t½ = 5730 years (NNDC)
        Isotope {
            nucleus: Nucleus::new(6, 14).unwrap_or_else(|_| Nucleus::carbon_12()),
            name: String::from("C-14"),
            half_life_seconds: 5730.0 * seconds_per_year,
            primary_decay: DecayMode::BetaMinus,
        },
        // Potassium-40: t½ = 1.248e9 years
        Isotope {
            nucleus: Nucleus::new(19, 40).unwrap_or_else(|_| Nucleus::hydrogen_1()),
            name: String::from("K-40"),
            half_life_seconds: 1.248e9 * seconds_per_year,
            primary_decay: DecayMode::BetaMinus,
        },
        // Cobalt-60: t½ = 5.2714 years
        Isotope {
            nucleus: Nucleus::new(27, 60).unwrap_or_else(|_| Nucleus::iron_56()),
            name: String::from("Co-60"),
            half_life_seconds: 5.2714 * seconds_per_year,
            primary_decay: DecayMode::BetaMinus,
        },
        // Strontium-90: t½ = 28.79 years
        Isotope {
            nucleus: Nucleus::new(38, 90).unwrap_or_else(|_| Nucleus::iron_56()),
            name: String::from("Sr-90"),
            half_life_seconds: 28.79 * seconds_per_year,
            primary_decay: DecayMode::BetaMinus,
        },
        // Iodine-131: t½ = 8.0197 days
        Isotope {
            nucleus: Nucleus::new(53, 131).unwrap_or_else(|_| Nucleus::iron_56()),
            name: String::from("I-131"),
            half_life_seconds: 8.0252 * seconds_per_day,
            primary_decay: DecayMode::BetaMinus,
        },
        // Cesium-137: t½ = 30.17 years
        Isotope {
            nucleus: Nucleus::new(55, 137).unwrap_or_else(|_| Nucleus::iron_56()),
            name: String::from("Cs-137"),
            half_life_seconds: 30.05 * seconds_per_year,
            primary_decay: DecayMode::BetaMinus,
        },
        // Radon-222: t½ = 3.8235 days
        Isotope {
            nucleus: Nucleus::new(86, 222).unwrap_or_else(|_| Nucleus::uranium_238()),
            name: String::from("Rn-222"),
            half_life_seconds: 3.8222 * seconds_per_day,
            primary_decay: DecayMode::Alpha,
        },
        // Radium-226: t½ = 1600 years
        Isotope {
            nucleus: Nucleus::new(88, 226).unwrap_or_else(|_| Nucleus::uranium_238()),
            name: String::from("Ra-226"),
            half_life_seconds: 1600.0 * seconds_per_year,
            primary_decay: DecayMode::Alpha,
        },
        // Thorium-230: t½ = 75380 years
        Isotope {
            nucleus: Nucleus::new(90, 230).unwrap_or_else(|_| Nucleus::uranium_238()),
            name: String::from("Th-230"),
            half_life_seconds: 75_380.0 * seconds_per_year,
            primary_decay: DecayMode::Alpha,
        },
        // Thorium-234: t½ = 24.10 days
        Isotope {
            nucleus: Nucleus::new(90, 234).unwrap_or_else(|_| Nucleus::uranium_238()),
            name: String::from("Th-234"),
            half_life_seconds: 24.10 * seconds_per_day,
            primary_decay: DecayMode::BetaMinus,
        },
        // Protactinium-234m: t½ = 1.17 minutes
        Isotope {
            nucleus: Nucleus::new(91, 234).unwrap_or_else(|_| Nucleus::uranium_238()),
            name: String::from("Pa-234m"),
            half_life_seconds: 1.17 * seconds_per_minute,
            primary_decay: DecayMode::BetaMinus,
        },
        // Uranium-234: t½ = 245500 years
        Isotope {
            nucleus: Nucleus::new(92, 234).unwrap_or_else(|_| Nucleus::uranium_235()),
            name: String::from("U-234"),
            half_life_seconds: 245_250.0 * seconds_per_year,
            primary_decay: DecayMode::Alpha,
        },
        // Uranium-235: t½ = 7.04e8 years
        Isotope {
            nucleus: Nucleus::uranium_235(),
            name: String::from("U-235"),
            half_life_seconds: 7.04e8 * seconds_per_year,
            primary_decay: DecayMode::Alpha,
        },
        // Uranium-238: t½ = 4.468e9 years (NNDC)
        Isotope {
            nucleus: Nucleus::uranium_238(),
            name: String::from("U-238"),
            half_life_seconds: 4.468e9 * seconds_per_year,
            primary_decay: DecayMode::Alpha,
        },
        // Plutonium-239: t½ = 24110 years
        Isotope {
            nucleus: Nucleus::new(94, 239).unwrap_or_else(|_| Nucleus::uranium_238()),
            name: String::from("Pu-239"),
            half_life_seconds: 24_110.0 * seconds_per_year,
            primary_decay: DecayMode::Alpha,
        },
        // Americium-241: t½ = 432.2 years
        Isotope {
            nucleus: Nucleus::new(95, 241).unwrap_or_else(|_| Nucleus::uranium_238()),
            name: String::from("Am-241"),
            half_life_seconds: 432.6 * seconds_per_year,
            primary_decay: DecayMode::Alpha,
        },
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
}
