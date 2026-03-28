//! Radioactive decay: decay modes, half-lives, activity, decay chains.

extern crate alloc;
use alloc::string::String;
use alloc::vec::Vec;

use crate::error::TanmatraError;
use crate::nucleus::Nucleus;
use serde::{Deserialize, Serialize};

/// Modes of radioactive decay.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum DecayMode {
    /// Alpha decay: emits He-4 nucleus (Z-2, A-4).
    Alpha,
    /// Beta-minus decay: neutron -> proton + electron + antineutrino (Z+1, A).
    BetaMinus,
    /// Beta-plus decay: proton -> neutron + positron + neutrino (Z-1, A).
    BetaPlus,
    /// Gamma decay: emits photon, no change in Z or A.
    Gamma,
    /// Electron capture: proton + electron -> neutron + neutrino (Z-1, A).
    ElectronCapture,
    /// Spontaneous fission.
    Fission,
}

/// Computes the decay constant lambda from a half-life in seconds.
///
/// `lambda = ln(2) / t_half`
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidEnergy`] if half-life is not positive and finite.
pub fn decay_constant(half_life_s: f64) -> Result<f64, TanmatraError> {
    if half_life_s <= 0.0 || half_life_s.is_nan() || half_life_s.is_infinite() {
        return Err(TanmatraError::InvalidEnergy {
            reason: String::from("half-life must be positive and finite"),
        });
    }
    Ok(core::f64::consts::LN_2 / half_life_s)
}

/// Fraction of atoms remaining after time t given decay constant lambda.
///
/// `N(t)/N(0) = e^(-lambda * t)`
#[must_use]
#[inline]
pub fn remaining_fraction(lambda: f64, time_s: f64) -> f64 {
    libm::exp(-lambda * time_s)
}

/// Activity in becquerels (decays per second).
///
/// `A = lambda * N`
#[must_use]
#[inline]
pub fn activity_bq(n_atoms: f64, lambda: f64) -> f64 {
    lambda * n_atoms
}

/// Applies alpha decay to a nucleus: Z-2, A-4.
///
/// # Errors
///
/// Returns error if the resulting nucleus would be invalid.
pub fn alpha_decay(z: u16, a: u16) -> Result<Nucleus, TanmatraError> {
    if z < 2 || a < 4 {
        return Err(TanmatraError::DecayFailed {
            reason: String::from("nucleus too light for alpha decay"),
        });
    }
    let new_z = z - 2;
    let new_a = a - 4;
    if new_a < new_z {
        return Err(TanmatraError::DecayFailed {
            reason: String::from("alpha decay would produce invalid nucleus"),
        });
    }
    Ok(Nucleus { z: new_z, a: new_a })
}

/// Applies beta-minus decay: Z+1, A unchanged.
///
/// # Errors
///
/// Returns error if Z+1 > A.
pub fn beta_minus_decay(z: u16, a: u16) -> Result<Nucleus, TanmatraError> {
    let new_z = z + 1;
    if new_z > a {
        return Err(TanmatraError::DecayFailed {
            reason: String::from("beta-minus decay would produce Z > A"),
        });
    }
    Ok(Nucleus { z: new_z, a })
}

/// Applies beta-plus decay: Z-1, A unchanged.
///
/// # Errors
///
/// Returns error if Z is 0.
pub fn beta_plus_decay(z: u16, a: u16) -> Result<Nucleus, TanmatraError> {
    if z == 0 {
        return Err(TanmatraError::DecayFailed {
            reason: String::from("cannot beta-plus decay a nucleus with Z=0"),
        });
    }
    Ok(Nucleus { z: z - 1, a })
}

/// Seconds in a Julian year (365.25 days).
const SECONDS_PER_YEAR: f64 = 365.25 * 24.0 * 3600.0;
/// Seconds in a day.
const SECONDS_PER_DAY: f64 = 24.0 * 3600.0;

/// A known radioactive isotope with experimentally measured properties.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct KnownIsotope {
    /// The nucleus.
    pub nucleus: Nucleus,
    /// Half-life in seconds.
    pub half_life_s: f64,
    /// Primary decay mode.
    pub primary_decay: DecayMode,
}

/// Carbon-14: half-life = 5730 years, beta-minus.
#[must_use]
pub fn carbon14() -> KnownIsotope {
    KnownIsotope {
        nucleus: Nucleus { z: 6, a: 14 },
        half_life_s: 5730.0 * SECONDS_PER_YEAR,
        primary_decay: DecayMode::BetaMinus,
    }
}

/// Potassium-40: half-life = 1.248e9 years, beta-minus.
#[must_use]
pub fn potassium40() -> KnownIsotope {
    KnownIsotope {
        nucleus: Nucleus { z: 19, a: 40 },
        half_life_s: 1.248e9 * SECONDS_PER_YEAR,
        primary_decay: DecayMode::BetaMinus,
    }
}

/// Cobalt-60: half-life = 5.271 years, beta-minus.
#[must_use]
pub fn cobalt60() -> KnownIsotope {
    KnownIsotope {
        nucleus: Nucleus { z: 27, a: 60 },
        half_life_s: 5.271 * SECONDS_PER_YEAR,
        primary_decay: DecayMode::BetaMinus,
    }
}

/// Iodine-131: half-life = 8.025 days, beta-minus.
#[must_use]
pub fn iodine131() -> KnownIsotope {
    KnownIsotope {
        nucleus: Nucleus { z: 53, a: 131 },
        half_life_s: 8.025 * SECONDS_PER_DAY,
        primary_decay: DecayMode::BetaMinus,
    }
}

/// Cesium-137: half-life = 30.17 years, beta-minus.
#[must_use]
pub fn cesium137() -> KnownIsotope {
    KnownIsotope {
        nucleus: Nucleus { z: 55, a: 137 },
        half_life_s: 30.17 * SECONDS_PER_YEAR,
        primary_decay: DecayMode::BetaMinus,
    }
}

/// Radium-226: half-life = 1600 years, alpha.
#[must_use]
pub fn radium226() -> KnownIsotope {
    KnownIsotope {
        nucleus: Nucleus { z: 88, a: 226 },
        half_life_s: 1600.0 * SECONDS_PER_YEAR,
        primary_decay: DecayMode::Alpha,
    }
}

/// Uranium-235: half-life = 7.04e8 years, alpha.
#[must_use]
pub fn uranium235() -> KnownIsotope {
    KnownIsotope {
        nucleus: Nucleus { z: 92, a: 235 },
        half_life_s: 7.04e8 * SECONDS_PER_YEAR,
        primary_decay: DecayMode::Alpha,
    }
}

/// Uranium-238: half-life = 4.468e9 years, alpha.
#[must_use]
pub fn uranium238() -> KnownIsotope {
    KnownIsotope {
        nucleus: Nucleus { z: 92, a: 238 },
        half_life_s: 4.468e9 * SECONDS_PER_YEAR,
        primary_decay: DecayMode::Alpha,
    }
}

/// Plutonium-239: half-life = 24110 years, alpha.
#[must_use]
pub fn plutonium239() -> KnownIsotope {
    KnownIsotope {
        nucleus: Nucleus { z: 94, a: 239 },
        half_life_s: 24110.0 * SECONDS_PER_YEAR,
        primary_decay: DecayMode::Alpha,
    }
}

/// Polonium-210: half-life = 138.4 days, alpha.
#[must_use]
pub fn polonium210() -> KnownIsotope {
    KnownIsotope {
        nucleus: Nucleus { z: 84, a: 210 },
        half_life_s: 138.4 * SECONDS_PER_DAY,
        primary_decay: DecayMode::Alpha,
    }
}

/// Follows a decay chain for the given number of steps, applying the
/// primary decay mode at each step.
///
/// Returns the sequence of daughter nuclei (not including the parent).
///
/// # Errors
///
/// Returns error if any decay step fails or if the decay mode is fission.
pub fn decay_chain(isotope: &KnownIsotope, steps: usize) -> Result<Vec<Nucleus>, TanmatraError> {
    let mut chain = Vec::with_capacity(steps);
    let mut current_z = isotope.nucleus.z;
    let mut current_a = isotope.nucleus.a;
    let mode = isotope.primary_decay;

    for _ in 0..steps {
        let daughter = match mode {
            DecayMode::Alpha => alpha_decay(current_z, current_a)?,
            DecayMode::BetaMinus => beta_minus_decay(current_z, current_a)?,
            DecayMode::BetaPlus | DecayMode::ElectronCapture => {
                beta_plus_decay(current_z, current_a)?
            }
            DecayMode::Gamma => Nucleus {
                z: current_z,
                a: current_a,
            },
            DecayMode::Fission => {
                return Err(TanmatraError::DecayFailed {
                    reason: String::from("fission products are not deterministic"),
                });
            }
        };
        chain.push(daughter);
        current_z = daughter.z;
        current_a = daughter.a;
    }

    Ok(chain)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn c14_half_life() {
        let c14 = carbon14();
        let expected_s = 5730.0 * SECONDS_PER_YEAR;
        assert!((c14.half_life_s - expected_s).abs() < 1.0);
    }

    #[test]
    fn c14_decay_constant() {
        let c14 = carbon14();
        let lambda = decay_constant(c14.half_life_s).unwrap();
        // lambda ~ 3.836e-12 /s
        let expected = 3.836e-12;
        let rel_err = libm::fabs(lambda - expected) / expected;
        assert!(rel_err < 0.01, "lambda = {lambda}, expected ~{expected}");
    }

    #[test]
    fn remaining_fraction_one_half_life() {
        let lambda = decay_constant(100.0).unwrap();
        let frac = remaining_fraction(lambda, 100.0);
        assert!((frac - 0.5).abs() < 1e-10);
    }

    #[test]
    fn alpha_decay_u238() {
        let daughter = alpha_decay(92, 238).unwrap();
        assert_eq!(daughter.z, 90); // Thorium
        assert_eq!(daughter.a, 234);
    }

    #[test]
    fn beta_minus_c14() {
        let daughter = beta_minus_decay(6, 14).unwrap();
        assert_eq!(daughter.z, 7); // Nitrogen
        assert_eq!(daughter.a, 14);
    }

    #[test]
    fn beta_plus_decay_basic() {
        let daughter = beta_plus_decay(8, 15).unwrap();
        assert_eq!(daughter.z, 7); // Nitrogen-15
        assert_eq!(daughter.a, 15);
    }

    #[test]
    fn alpha_decay_too_light() {
        assert!(alpha_decay(1, 3).is_err());
    }

    #[test]
    fn decay_chain_u238_alpha() {
        let u238 = uranium238();
        let chain = decay_chain(&u238, 3).unwrap();
        assert_eq!(chain.len(), 3);
        // First alpha: U-238 -> Th-234
        assert_eq!(chain[0].z, 90);
        assert_eq!(chain[0].a, 234);
        // Second alpha: Th-234 -> Ra-230
        assert_eq!(chain[1].z, 88);
        assert_eq!(chain[1].a, 230);
        // Third alpha: Ra-230 -> Rn-226
        assert_eq!(chain[2].z, 86);
        assert_eq!(chain[2].a, 226);
    }

    #[test]
    fn serde_roundtrip_decay_mode() {
        let modes = [
            DecayMode::Alpha,
            DecayMode::BetaMinus,
            DecayMode::BetaPlus,
            DecayMode::Gamma,
            DecayMode::ElectronCapture,
            DecayMode::Fission,
        ];
        for mode in &modes {
            let json = serde_json::to_string(mode).unwrap();
            let back: DecayMode = serde_json::from_str(&json).unwrap();
            assert_eq!(*mode, back);
        }
    }

    #[test]
    fn serde_roundtrip_known_isotope() {
        let c14 = carbon14();
        let json = serde_json::to_string(&c14).unwrap();
        let back: KnownIsotope = serde_json::from_str(&json).unwrap();
        assert_eq!(c14.nucleus, back.nucleus);
        assert_eq!(c14.primary_decay, back.primary_decay);
    }

    #[test]
    fn activity_basic() {
        let lambda = 0.001;
        let n = 1e10;
        let a = activity_bq(n, lambda);
        assert!((a - 1e7).abs() < 1.0);
    }
}
