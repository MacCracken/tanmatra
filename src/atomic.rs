//! Electron structure: orbital types, quantum numbers, electron configuration,
//! spectral lines, and ionization energies.

extern crate alloc;
use alloc::string::String;
use alloc::vec::Vec;

use crate::constants::RYDBERG_CONST;
use crate::error::TanmatraError;
use serde::{Deserialize, Serialize};

/// Orbital angular momentum types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum OrbitalType {
    /// l=0, sharp.
    S,
    /// l=1, principal.
    P,
    /// l=2, diffuse.
    D,
    /// l=3, fundamental.
    F,
}

impl OrbitalType {
    /// Returns the angular momentum quantum number l for this orbital.
    #[must_use]
    #[inline]
    pub fn l(self) -> u8 {
        match self {
            Self::S => 0,
            Self::P => 1,
            Self::D => 2,
            Self::F => 3,
        }
    }

    /// Creates an `OrbitalType` from the angular momentum quantum number.
    ///
    /// # Errors
    ///
    /// Returns [`TanmatraError::ComputationError`] if l is not 0-3.
    pub fn from_l(l: u8) -> Result<Self, TanmatraError> {
        match l {
            0 => Ok(Self::S),
            1 => Ok(Self::P),
            2 => Ok(Self::D),
            3 => Ok(Self::F),
            _ => Err(TanmatraError::ComputationError {
                reason: String::from("l must be 0-3 (s, p, d, f)"),
            }),
        }
    }
}

/// A complete set of quantum numbers for an electron.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct QuantumNumbers {
    /// Principal quantum number (n >= 1).
    pub n: u8,
    /// Angular momentum quantum number (0 <= l < n).
    pub l: u8,
    /// Magnetic quantum number (-l <= ml <= l).
    pub ml: i8,
    /// Spin quantum number (+1/2 or -1/2).
    pub ms: f64,
}

impl QuantumNumbers {
    /// Creates and validates a set of quantum numbers.
    ///
    /// # Errors
    ///
    /// Returns error if the quantum numbers violate selection rules.
    pub fn new(n: u8, l: u8, ml: i8, ms: f64) -> Result<Self, TanmatraError> {
        if n == 0 {
            return Err(TanmatraError::ComputationError {
                reason: String::from("n must be >= 1"),
            });
        }
        if l >= n {
            return Err(TanmatraError::ComputationError {
                reason: String::from("l must be < n"),
            });
        }
        let l_signed = l.cast_signed();
        if ml < -l_signed || ml > l_signed {
            return Err(TanmatraError::ComputationError {
                reason: String::from("ml must be in range [-l, l]"),
            });
        }
        if (ms - 0.5).abs() > f64::EPSILON && (ms + 0.5).abs() > f64::EPSILON {
            return Err(TanmatraError::ComputationError {
                reason: String::from("ms must be +0.5 or -0.5"),
            });
        }
        Ok(Self { n, l, ml, ms })
    }
}

/// Maximum electrons in a subshell with angular momentum l: 2*(2l+1).
#[must_use]
#[inline]
pub fn max_electrons_in_subshell(l: u8) -> u8 {
    2 * (2 * l + 1)
}

/// Maximum electrons in a shell with principal quantum number n: 2*n^2.
#[must_use]
#[inline]
pub fn max_electrons_in_shell(n: u8) -> u16 {
    2 * u16::from(n) * u16::from(n)
}

/// Aufbau filling order using the Madelung (n+l) rule.
/// Returns (n, l) pairs in filling order.
fn aufbau_order() -> Vec<(u8, u8)> {
    let mut subshells = Vec::new();
    for n in 1u8..=7 {
        for l in 0..n.min(4) {
            subshells.push((n, l));
        }
    }
    subshells.sort_by_key(|&(n, l)| (n + l, n));
    subshells
}

/// Computes the ground-state electron configuration for element with atomic number Z.
///
/// Returns a list of (n, `OrbitalType`, count) tuples.
///
/// Handles the well-known exceptions:
/// - Cr (Z=24): \[Ar\] 3d5 4s1 (half-filled d shell)
/// - Cu (Z=29): \[Ar\] 3d10 4s1 (fully-filled d shell)
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidElement`] if Z is 0 or > 118.
pub fn electron_configuration(z: u16) -> Result<Vec<(u8, OrbitalType, u8)>, TanmatraError> {
    if z == 0 || z > 118 {
        return Err(TanmatraError::InvalidElement { z });
    }

    let order = aufbau_order();
    let mut config: Vec<(u8, u8, u8)> = Vec::new(); // (n, l, count)
    #[allow(clippy::cast_possible_truncation)]
    let mut remaining = z as u8; // z <= 118, safe

    for &(n, l) in &order {
        if remaining == 0 {
            break;
        }
        let max_e = max_electrons_in_subshell(l);
        let fill = remaining.min(max_e);
        config.push((n, l, fill));
        remaining -= fill;
    }

    // Handle Cr (Z=24) and Cu (Z=29) exceptions.
    if z == 24 || z == 29 {
        let mut d_index = None;
        let mut s_index = None;
        for (i, &(n, l, _)) in config.iter().enumerate() {
            if n == 3 && l == 2 {
                d_index = Some(i);
            }
            if n == 4 && l == 0 {
                s_index = Some(i);
            }
        }
        if let (Some(di), Some(si)) = (d_index, s_index)
            && config[si].2 >= 1
        {
            config[si].2 -= 1;
            config[di].2 += 1;
        }
    }

    config
        .into_iter()
        .filter(|&(_, _, c)| c > 0)
        .map(|(n, l, c)| OrbitalType::from_l(l).map(|ot| (n, ot, c)))
        .collect()
}

/// Computes the wavelength (in nm) of a spectral line for a hydrogen-like
/// atom using the Rydberg formula.
///
/// `1/lambda = R_inf * Z^2 * (1/n_lower^2 - 1/n_upper^2)`
///
/// # Errors
///
/// Returns error if `n_upper` <= `n_lower`, or either is zero.
pub fn spectral_line_nm(z: u16, n_upper: u16, n_lower: u16) -> Result<f64, TanmatraError> {
    if n_lower == 0 || n_upper == 0 {
        return Err(TanmatraError::ComputationError {
            reason: String::from("quantum numbers must be > 0"),
        });
    }
    if n_upper <= n_lower {
        return Err(TanmatraError::ComputationError {
            reason: String::from("n_upper must be greater than n_lower"),
        });
    }

    let zf = f64::from(z);
    let nl = f64::from(n_lower);
    let nu = f64::from(n_upper);

    let inv_lambda = RYDBERG_CONST * zf * zf * (1.0 / (nl * nl) - 1.0 / (nu * nu));

    if inv_lambda <= 0.0 {
        return Err(TanmatraError::ComputationError {
            reason: String::from("computed inverse wavelength is non-positive"),
        });
    }

    Ok(1.0e9 / inv_lambda)
}

/// Computes visible Balmer series lines for hydrogen (transitions to n=2).
///
/// Returns the wavelength in nm for a transition from `n_upper` to n=2.
///
/// # Errors
///
/// Returns error if `n_upper` <= 2.
#[inline]
pub fn balmer_series(n_upper: u16) -> Result<f64, TanmatraError> {
    spectral_line_nm(1, n_upper, 2)
}

/// First ionization energy (in eV) for elements Z=1 to Z=36.
///
/// Values from NIST Atomic Spectra Database.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidElement`] if Z is out of range 1-36.
pub fn ionization_energy_ev(z: u16) -> Result<f64, TanmatraError> {
    #[rustfmt::skip]
    let table: [f64; 36] = [
        13.598, 24.587,  5.392,  9.323,  8.298, 11.260, 14.534, 13.618, // H-O
        17.423, 21.565,  5.139,  7.646,  5.986,  8.152, 10.487, 10.360, // F-S
        12.968, 15.760,  4.341,  6.113,  6.562,  6.828,  6.746,  6.767, // Cl-Cr
         7.434,  7.902,  7.881,  7.640,  7.726,  9.394,  5.999,  7.900, // Mn-Ge
         9.815,  9.752, 11.814, 14.000,                                   // As-Kr
    ];

    if z == 0 || z > 36 {
        return Err(TanmatraError::InvalidElement { z });
    }

    Ok(table[(z - 1) as usize])
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    #[test]
    fn orbital_type_l_values() {
        assert_eq!(OrbitalType::S.l(), 0);
        assert_eq!(OrbitalType::P.l(), 1);
        assert_eq!(OrbitalType::D.l(), 2);
        assert_eq!(OrbitalType::F.l(), 3);
    }

    #[test]
    fn max_electrons_subshell() {
        assert_eq!(max_electrons_in_subshell(0), 2);
        assert_eq!(max_electrons_in_subshell(1), 6);
        assert_eq!(max_electrons_in_subshell(2), 10);
        assert_eq!(max_electrons_in_subshell(3), 14);
    }

    #[test]
    fn max_electrons_shell() {
        assert_eq!(max_electrons_in_shell(1), 2);
        assert_eq!(max_electrons_in_shell(2), 8);
        assert_eq!(max_electrons_in_shell(3), 18);
        assert_eq!(max_electrons_in_shell(4), 32);
    }

    #[test]
    fn quantum_numbers_validation() {
        assert!(QuantumNumbers::new(1, 0, 0, 0.5).is_ok());
        assert!(QuantumNumbers::new(2, 1, -1, -0.5).is_ok());
        assert!(QuantumNumbers::new(0, 0, 0, 0.5).is_err());
        assert!(QuantumNumbers::new(1, 1, 0, 0.5).is_err());
        assert!(QuantumNumbers::new(2, 1, 2, 0.5).is_err());
        assert!(QuantumNumbers::new(1, 0, 0, 0.3).is_err());
    }

    #[test]
    fn hydrogen_config() {
        let config = electron_configuration(1).unwrap();
        assert_eq!(config, vec![(1, OrbitalType::S, 1)]);
    }

    #[test]
    fn helium_config() {
        let config = electron_configuration(2).unwrap();
        assert_eq!(config, vec![(1, OrbitalType::S, 2)]);
    }

    #[test]
    fn iron_config() {
        let config = electron_configuration(26).unwrap();
        let d_shell = config
            .iter()
            .find(|&&(n, ot, _)| n == 3 && ot == OrbitalType::D);
        let s_shell = config
            .iter()
            .find(|&&(n, ot, _)| n == 4 && ot == OrbitalType::S);
        assert_eq!(d_shell, Some(&(3, OrbitalType::D, 6)));
        assert_eq!(s_shell, Some(&(4, OrbitalType::S, 2)));
    }

    #[test]
    fn chromium_exception() {
        let config = electron_configuration(24).unwrap();
        let d_shell = config
            .iter()
            .find(|&&(n, ot, _)| n == 3 && ot == OrbitalType::D);
        let s_shell = config
            .iter()
            .find(|&&(n, ot, _)| n == 4 && ot == OrbitalType::S);
        assert_eq!(d_shell, Some(&(3, OrbitalType::D, 5)));
        assert_eq!(s_shell, Some(&(4, OrbitalType::S, 1)));
    }

    #[test]
    fn copper_exception() {
        let config = electron_configuration(29).unwrap();
        let d_shell = config
            .iter()
            .find(|&&(n, ot, _)| n == 3 && ot == OrbitalType::D);
        let s_shell = config
            .iter()
            .find(|&&(n, ot, _)| n == 4 && ot == OrbitalType::S);
        assert_eq!(d_shell, Some(&(3, OrbitalType::D, 10)));
        assert_eq!(s_shell, Some(&(4, OrbitalType::S, 1)));
    }

    #[test]
    fn hydrogen_balmer_h_alpha() {
        let lambda = balmer_series(3).unwrap();
        assert!(
            (lambda - 656.3).abs() < 0.5,
            "H-alpha = {lambda} nm, expected ~656.3"
        );
    }

    #[test]
    fn spectral_line_invalid() {
        assert!(spectral_line_nm(1, 2, 2).is_err());
        assert!(spectral_line_nm(1, 1, 2).is_err());
        assert!(spectral_line_nm(1, 0, 1).is_err());
    }

    #[test]
    fn ionization_energy_hydrogen() {
        let ie = ionization_energy_ev(1).unwrap();
        assert!((ie - 13.598).abs() < 0.001);
    }

    #[test]
    fn ionization_energy_noble_gases_high() {
        let he = ionization_energy_ev(2).unwrap();
        let ne = ionization_energy_ev(10).unwrap();
        let ar = ionization_energy_ev(18).unwrap();
        let kr = ionization_energy_ev(36).unwrap();

        assert!(he > ionization_energy_ev(3).unwrap());
        assert!(ne > ionization_energy_ev(11).unwrap());
        assert!(ar > ionization_energy_ev(19).unwrap());
        assert!(kr > ionization_energy_ev(31).unwrap());
    }

    #[test]
    fn ionization_energy_out_of_range() {
        assert!(ionization_energy_ev(0).is_err());
        assert!(ionization_energy_ev(37).is_err());
    }

    #[test]
    fn total_electrons_in_config() {
        for z in 1u16..=36 {
            let config = electron_configuration(z).unwrap();
            let total: u16 = config.iter().map(|&(_, _, c)| u16::from(c)).sum();
            assert_eq!(total, z, "Z={z} config has {total} electrons");
        }
    }

    #[test]
    fn serde_roundtrip_quantum_numbers() {
        let qn = QuantumNumbers::new(3, 2, -1, 0.5).unwrap();
        let json = serde_json::to_string(&qn).unwrap();
        let back: QuantumNumbers = serde_json::from_str(&json).unwrap();
        assert_eq!(qn, back);
    }

    #[test]
    fn serde_roundtrip_orbital_type() {
        let ot = OrbitalType::D;
        let json = serde_json::to_string(&ot).unwrap();
        let back: OrbitalType = serde_json::from_str(&json).unwrap();
        assert_eq!(ot, back);
    }
}
