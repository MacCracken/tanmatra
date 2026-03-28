//! Nuclear structure: binding energy, mass defect, nuclear radius.
//!
//! Uses the Bethe-Weizsacker semi-empirical mass formula (liquid drop model)
//! with standard coefficients.

use crate::constants::R0_FM;
use crate::error::TanmatraError;
use serde::{Deserialize, Serialize};

/// A nucleus characterized by atomic number Z and mass number A.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Nucleus {
    /// Atomic number (number of protons).
    pub z: u16,
    /// Mass number (number of protons + neutrons).
    pub a: u16,
}

impl Nucleus {
    /// Creates a new nucleus, validating that A >= Z and both are > 0.
    ///
    /// # Errors
    ///
    /// Returns [`TanmatraError::InvalidIsotope`] if A < Z or either is zero.
    pub fn new(z: u16, a: u16) -> Result<Self, TanmatraError> {
        if z == 0 || a == 0 || a < z {
            return Err(TanmatraError::InvalidIsotope { z, a });
        }
        Ok(Self { z, a })
    }

    /// Number of neutrons (N = A - Z).
    #[must_use]
    #[inline]
    pub fn n(&self) -> u16 {
        self.a - self.z
    }
}

// Bethe-Weizsacker coefficients (all in MeV).
const A_V: f64 = 15.67;
const A_S: f64 = 17.23;
const A_C: f64 = 0.714;
const A_A: f64 = 93.15;
const DELTA_0: f64 = 12.0;

/// Computes the nuclear binding energy (in `MeV`) using the Bethe-Weizsacker
/// semi-empirical mass formula.
///
/// `B = a_v*A - a_s*A^(2/3) - a_c*Z*(Z-1)/A^(1/3) - a_a*(A-2Z)^2/(4A) + delta`
///
/// where delta = +12/sqrt(A) for even-even, -12/sqrt(A) for odd-odd, 0 otherwise.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidIsotope`] if Z or A is zero, or A < Z.
pub fn binding_energy_mev(z: u16, a: u16) -> Result<f64, TanmatraError> {
    if z == 0 || a == 0 || a < z {
        return Err(TanmatraError::InvalidIsotope { z, a });
    }

    let af = f64::from(a);
    let zf = f64::from(z);

    let volume = A_V * af;
    let surface = A_S * libm::pow(af, 2.0 / 3.0);
    let coulomb = A_C * zf * (zf - 1.0) / libm::pow(af, 1.0 / 3.0);
    let asymmetry = A_A * (af - 2.0 * zf) * (af - 2.0 * zf) / (4.0 * af);

    let z_even = z.is_multiple_of(2);
    let n_even = (a - z).is_multiple_of(2);
    let delta = if z_even && n_even {
        DELTA_0 / libm::sqrt(af)
    } else if !z_even && !n_even {
        -DELTA_0 / libm::sqrt(af)
    } else {
        0.0
    };

    let b = volume - surface - coulomb - asymmetry + delta;
    Ok(b)
}

/// Binding energy per nucleon (in `MeV`).
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidIsotope`] if Z or A is invalid.
pub fn binding_energy_per_nucleon(z: u16, a: u16) -> Result<f64, TanmatraError> {
    let b = binding_energy_mev(z, a)?;
    Ok(b / f64::from(a))
}

/// Mass defect (in `MeV`): equals the binding energy in the semi-empirical model.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidIsotope`] if Z or A is invalid.
pub fn mass_defect_mev(z: u16, a: u16) -> Result<f64, TanmatraError> {
    binding_energy_mev(z, a)
}

/// Nuclear radius (in fm) using the empirical formula `R = R0 * A^(1/3)`.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidIsotope`] if A is zero.
#[inline]
pub fn nuclear_radius_fm(a: u16) -> Result<f64, TanmatraError> {
    if a == 0 {
        return Err(TanmatraError::InvalidIsotope { z: 0, a: 0 });
    }
    Ok(R0_FM * libm::cbrt(f64::from(a)))
}

/// Checks whether a number is a nuclear magic number.
///
/// Magic numbers correspond to completely filled nuclear shells and confer
/// extra stability: 2, 8, 20, 28, 50, 82, 126.
#[must_use]
#[inline]
pub fn is_magic_number(n: u16) -> bool {
    matches!(n, 2 | 8 | 20 | 28 | 50 | 82 | 126)
}

// Preset nuclei.

/// Hydrogen-1 (protium).
pub const HYDROGEN: Nucleus = Nucleus { z: 1, a: 1 };
/// Helium-4 (alpha particle).
pub const HELIUM4: Nucleus = Nucleus { z: 2, a: 4 };
/// Carbon-12.
pub const CARBON12: Nucleus = Nucleus { z: 6, a: 12 };
/// Oxygen-16.
pub const OXYGEN16: Nucleus = Nucleus { z: 8, a: 16 };
/// Iron-56 (peak of binding energy curve).
pub const IRON56: Nucleus = Nucleus { z: 26, a: 56 };
/// Uranium-235 (fissile).
pub const URANIUM235: Nucleus = Nucleus { z: 92, a: 235 };
/// Uranium-238 (most common uranium isotope).
pub const URANIUM238: Nucleus = Nucleus { z: 92, a: 238 };

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nucleus_neutron_count() {
        let u238 = Nucleus { z: 92, a: 238 };
        assert_eq!(u238.n(), 146);
    }

    #[test]
    fn invalid_nucleus_a_less_than_z() {
        assert!(Nucleus::new(10, 5).is_err());
    }

    #[test]
    fn invalid_nucleus_zero() {
        assert!(Nucleus::new(0, 0).is_err());
    }

    #[test]
    fn fe56_binding_energy_per_nucleon() {
        let bea = binding_energy_per_nucleon(26, 56).unwrap();
        assert!(bea > 8.5 && bea < 9.1, "Fe-56 B/A = {bea}");
    }

    #[test]
    fn he4_binding_energy() {
        let b = binding_energy_mev(2, 4).unwrap();
        assert!(b > 24.0 && b < 32.0, "He-4 B = {b}");
    }

    #[test]
    fn nuclear_radius_he4() {
        let r = nuclear_radius_fm(4).unwrap();
        assert!((r - 1.905).abs() < 0.01, "He-4 R = {r}");
    }

    #[test]
    fn magic_numbers() {
        let magic = [2, 8, 20, 28, 50, 82, 126];
        for m in &magic {
            assert!(is_magic_number(*m), "{m} should be magic");
        }
        assert!(!is_magic_number(3));
        assert!(!is_magic_number(100));
    }

    #[test]
    fn binding_energy_increases_with_a() {
        let b12 = binding_energy_mev(6, 12).unwrap();
        let b56 = binding_energy_mev(26, 56).unwrap();
        assert!(b56 > b12);
    }

    #[test]
    fn serde_roundtrip_nucleus() {
        let n = Nucleus { z: 26, a: 56 };
        let json = serde_json::to_string(&n).unwrap();
        let back: Nucleus = serde_json::from_str(&json).unwrap();
        assert_eq!(n, back);
    }
}
