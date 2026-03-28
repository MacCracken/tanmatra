//! Atomic structure: quantum numbers, electron configurations, and spectral lines.
//!
//! Implements the hydrogen-like spectral line formula (Rydberg), electron
//! configurations via the Aufbau principle with Madelung's rule, and
//! ionization energies from NIST data.

use crate::constants::RYDBERG;
use crate::error::TanmatraError;
use alloc::string::String;
use alloc::vec::Vec;
use serde::{Deserialize, Serialize};

/// Orbital angular momentum quantum number types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum OrbitalType {
    /// l = 0
    S,
    /// l = 1
    P,
    /// l = 2
    D,
    /// l = 3
    F,
}

impl OrbitalType {
    /// Returns the angular momentum quantum number l.
    #[must_use]
    pub const fn l(self) -> u32 {
        match self {
            Self::S => 0,
            Self::P => 1,
            Self::D => 2,
            Self::F => 3,
        }
    }

    /// Returns the maximum number of electrons in this orbital type.
    #[must_use]
    pub const fn max_electrons(self) -> u32 {
        2 * (2 * self.l() + 1)
    }

    /// Returns the letter symbol.
    #[must_use]
    pub const fn symbol(self) -> char {
        match self {
            Self::S => 's',
            Self::P => 'p',
            Self::D => 'd',
            Self::F => 'f',
        }
    }
}

/// A set of quantum numbers (n, l, ml, ms) describing an electron state.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct QuantumNumbers {
    /// Principal quantum number (n >= 1).
    pub n: u32,
    /// Orbital angular momentum quantum number (0 <= l < n).
    pub l: u32,
    /// Magnetic quantum number (-l <= ml <= l).
    pub ml: i32,
    /// Spin magnetic quantum number (+1/2 or -1/2, stored as +1 or -1).
    pub ms: i32,
}

impl QuantumNumbers {
    /// Creates and validates a set of quantum numbers.
    ///
    /// # Errors
    ///
    /// Returns [`TanmatraError::InvalidQuantumNumbers`] if any constraint is violated.
    pub fn new(n: u32, l: u32, ml: i32, ms: i32) -> Result<Self, TanmatraError> {
        if n == 0 {
            return Err(TanmatraError::InvalidQuantumNumbers(String::from(
                "n must be >= 1",
            )));
        }
        if l >= n {
            return Err(TanmatraError::InvalidQuantumNumbers(alloc::format!(
                "l={l} must be < n={n}"
            )));
        }
        if ml.unsigned_abs() > l {
            return Err(TanmatraError::InvalidQuantumNumbers(alloc::format!(
                "ml={ml} must satisfy |ml| <= l={l}"
            )));
        }
        if ms != 1 && ms != -1 {
            return Err(TanmatraError::InvalidQuantumNumbers(alloc::format!(
                "ms={ms} must be +1 or -1 (representing +1/2 or -1/2)"
            )));
        }
        Ok(Self { n, l, ml, ms })
    }
}

/// An orbital filling entry for electron configuration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct OrbitalFilling {
    /// Principal quantum number.
    pub n: u32,
    /// Orbital type.
    pub orbital: OrbitalType,
    /// Number of electrons in this subshell.
    pub electrons: u32,
}

/// Aufbau filling order using Madelung's rule (n+l, then n).
///
/// Returns subshells in filling order: 1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, ...
const FILLING_ORDER: [(u32, OrbitalType); 20] = [
    (1, OrbitalType::S),
    (2, OrbitalType::S),
    (2, OrbitalType::P),
    (3, OrbitalType::S),
    (3, OrbitalType::P),
    (4, OrbitalType::S),
    (3, OrbitalType::D),
    (4, OrbitalType::P),
    (5, OrbitalType::S),
    (4, OrbitalType::D),
    (5, OrbitalType::P),
    (6, OrbitalType::S),
    (4, OrbitalType::F),
    (5, OrbitalType::D),
    (6, OrbitalType::P),
    (7, OrbitalType::S),
    (5, OrbitalType::F),
    (6, OrbitalType::D),
    (7, OrbitalType::P),
    (8, OrbitalType::S),
];

/// Returns the electron configuration for element with atomic number Z.
///
/// Follows the Aufbau principle with Madelung's rule, including the
/// well-known exceptions for chromium (Z=24) and copper (Z=29).
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidAtomicNumber`] if `z` is 0 or > 118.
pub fn electron_configuration(z: u32) -> Result<Vec<OrbitalFilling>, TanmatraError> {
    if z == 0 || z > 118 {
        return Err(TanmatraError::InvalidAtomicNumber(z));
    }

    let mut config = Vec::new();
    let mut remaining = z;

    for &(n, orbital) in &FILLING_ORDER {
        if remaining == 0 {
            break;
        }
        let max = orbital.max_electrons();
        let fill = if remaining >= max { max } else { remaining };
        config.push(OrbitalFilling {
            n,
            orbital,
            electrons: fill,
        });
        remaining -= fill;
    }

    // Apply known exceptions
    apply_exceptions(z, &mut config);

    Ok(config)
}

/// Applies known electron configuration exceptions.
///
/// Chromium (Z=24): [Ar] 3d5 4s1 instead of [Ar] 3d4 4s2
/// Copper (Z=29): [Ar] 3d10 4s1 instead of [Ar] 3d9 4s2
/// Molybdenum (Z=42): [Kr] 4d5 5s1 instead of [Kr] 4d4 5s2
/// Silver (Z=47): [Kr] 4d10 5s1 instead of [Kr] 4d9 5s2
/// Gold (Z=79): [Xe] 4f14 5d10 6s1 instead of [Xe] 4f14 5d9 6s2
fn apply_exceptions(z: u32, config: &mut [OrbitalFilling]) {
    match z {
        24 => {
            // Cr: move one electron from 4s to 3d
            for entry in config.iter_mut() {
                if entry.n == 4 && entry.orbital == OrbitalType::S {
                    entry.electrons = 1;
                }
                if entry.n == 3 && entry.orbital == OrbitalType::D {
                    entry.electrons = 5;
                }
            }
        }
        29 => {
            // Cu: move one electron from 4s to 3d
            for entry in config.iter_mut() {
                if entry.n == 4 && entry.orbital == OrbitalType::S {
                    entry.electrons = 1;
                }
                if entry.n == 3 && entry.orbital == OrbitalType::D {
                    entry.electrons = 10;
                }
            }
        }
        42 => {
            // Mo: move one electron from 5s to 4d
            for entry in config.iter_mut() {
                if entry.n == 5 && entry.orbital == OrbitalType::S {
                    entry.electrons = 1;
                }
                if entry.n == 4 && entry.orbital == OrbitalType::D {
                    entry.electrons = 5;
                }
            }
        }
        47 => {
            // Ag: move one electron from 5s to 4d
            for entry in config.iter_mut() {
                if entry.n == 5 && entry.orbital == OrbitalType::S {
                    entry.electrons = 1;
                }
                if entry.n == 4 && entry.orbital == OrbitalType::D {
                    entry.electrons = 10;
                }
            }
        }
        79 => {
            // Au: move one electron from 6s to 5d
            for entry in config.iter_mut() {
                if entry.n == 6 && entry.orbital == OrbitalType::S {
                    entry.electrons = 1;
                }
                if entry.n == 5 && entry.orbital == OrbitalType::D {
                    entry.electrons = 10;
                }
            }
        }
        _ => {}
    }
}

/// Formats an electron configuration as a string.
///
/// Example: `"1s2 2s2 2p6 3s2 3p6 4s2 3d6"` for iron (Z=26).
#[must_use]
pub fn format_configuration(config: &[OrbitalFilling]) -> String {
    let mut parts = Vec::new();
    for entry in config {
        parts.push(alloc::format!(
            "{}{}{}",
            entry.n,
            entry.orbital.symbol(),
            entry.electrons
        ));
    }
    let mut result = String::new();
    for (i, part) in parts.iter().enumerate() {
        if i > 0 {
            result.push(' ');
        }
        result.push_str(part);
    }
    result
}

/// Formats an electron configuration with noble gas core notation.
///
/// Example: `"[Ar] 3d6 4s2"` for iron (Z=26).
#[must_use]
pub fn format_configuration_short(config: &[OrbitalFilling], z: u32) -> String {
    // Noble gas cores: He=2, Ne=10, Ar=18, Kr=36, Xe=54, Rn=86
    let (core_symbol, core_z) = if z > 86 {
        ("[Rn]", 86u32)
    } else if z > 54 {
        ("[Xe]", 54)
    } else if z > 36 {
        ("[Kr]", 36)
    } else if z > 18 {
        ("[Ar]", 18)
    } else if z > 10 {
        ("[Ne]", 10)
    } else if z > 2 {
        ("[He]", 2)
    } else {
        return format_configuration(config);
    };

    // Count electrons in the core
    let mut core_electrons = 0u32;
    let mut valence_start = 0;
    for (i, entry) in config.iter().enumerate() {
        core_electrons += entry.electrons;
        if core_electrons >= core_z {
            valence_start = i + 1;
            break;
        }
    }

    if valence_start >= config.len() {
        return String::from(core_symbol);
    }

    let valence = &config[valence_start..];
    let valence_str = format_configuration(valence);

    alloc::format!("{core_symbol} {valence_str}")
}

/// Calculates the wavelength of a spectral line in nanometers using
/// the Rydberg formula for hydrogen-like atoms.
///
/// 1/lambda = R_inf * Z^2 * |1/n1^2 - 1/n2^2|
///
/// where n1 < n2 (n1 is the lower energy level).
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] if n1 or n2 is 0, or n1 == n2.
#[inline]
pub fn spectral_line_nm(z: u32, n1: u32, n2: u32) -> Result<f64, TanmatraError> {
    if n1 == 0 || n2 == 0 {
        return Err(TanmatraError::InvalidQuantumNumbers(String::from(
            "quantum numbers must be >= 1",
        )));
    }
    if n1 == n2 {
        return Err(TanmatraError::InvalidQuantumNumbers(String::from(
            "n1 and n2 must differ for a transition",
        )));
    }

    let (lower, upper) = if n1 < n2 { (n1, n2) } else { (n2, n1) };

    let z_f = z as f64;
    let inv_lambda = RYDBERG
        * z_f
        * z_f
        * (1.0 / (lower as f64 * lower as f64) - 1.0 / (upper as f64 * upper as f64));

    // Convert from m^-1 to nm
    Ok(1.0e9 / inv_lambda)
}

/// NIST ionization energies for Z=1 to Z=36 in eV.
///
/// Source: NIST Atomic Spectra Database, Ionization Energies.
/// <https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html>
const IONIZATION_ENERGIES_EV: [f64; 36] = [
    13.598_44, // H  (Z=1)
    24.587_4,  // He (Z=2)
    5.391_71,  // Li (Z=3)
    9.322_7,   // Be (Z=4)
    8.298_0,   // B  (Z=5)
    11.260_3,  // C  (Z=6)
    14.534_1,  // N  (Z=7)
    13.618_1,  // O  (Z=8)
    17.422_8,  // F  (Z=9)
    21.564_5,  // Ne (Z=10)
    5.139_08,  // Na (Z=11)
    7.646_2,   // Mg (Z=12)
    5.985_77,  // Al (Z=13)
    8.151_7,   // Si (Z=14)
    10.486_7,  // P  (Z=15)
    10.360_0,  // S  (Z=16)
    12.967_6,  // Cl (Z=17)
    15.759_6,  // Ar (Z=18)
    4.340_66,  // K  (Z=19)
    6.113_2,   // Ca (Z=20)
    6.561_5,   // Sc (Z=21)
    6.828_1,   // Ti (Z=22)
    6.746_2,   // V  (Z=23)
    6.766_5,   // Cr (Z=24)
    7.434_0,   // Mn (Z=25)
    7.902_4,   // Fe (Z=26)
    7.881_0,   // Co (Z=27)
    7.639_8,   // Ni (Z=28)
    7.726_4,   // Cu (Z=29)
    9.394_2,   // Zn (Z=30)
    5.999_3,   // Ga (Z=31)
    7.899_4,   // Ge (Z=32)
    9.788_6,   // As (Z=33)
    9.752_4,   // Se (Z=34)
    11.813_8,  // Br (Z=35)
    13.999_6,  // Kr (Z=36)
];

/// Returns the first ionization energy in eV for the given atomic number.
///
/// Covers Z=1 (hydrogen) through Z=36 (krypton) using NIST values.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidAtomicNumber`] if Z is 0 or > 36.
#[inline]
pub fn ionization_energy_ev(z: u32) -> Result<f64, TanmatraError> {
    if z == 0 || z > 36 {
        return Err(TanmatraError::InvalidAtomicNumber(z));
    }
    Ok(IONIZATION_ENERGIES_EV[(z - 1) as usize])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn h_alpha_wavelength() {
        // H-alpha: n=3 -> n=2, hydrogen (Z=1)
        // Expected: ~656.3 nm
        let lambda = spectral_line_nm(1, 2, 3).unwrap();
        assert!(
            (lambda - 656.3).abs() < 1.0,
            "H-alpha wavelength={lambda} nm"
        );
    }

    #[test]
    fn lyman_alpha() {
        // Lyman-alpha: n=2 -> n=1, hydrogen
        // Expected: ~121.6 nm
        let lambda = spectral_line_nm(1, 1, 2).unwrap();
        assert!(
            (lambda - 121.6).abs() < 0.5,
            "Lyman-alpha wavelength={lambda} nm"
        );
    }

    #[test]
    fn spectral_line_invalid() {
        assert!(spectral_line_nm(1, 0, 2).is_err());
        assert!(spectral_line_nm(1, 2, 2).is_err());
    }

    #[test]
    fn iron_electron_config() {
        let config = electron_configuration(26).unwrap();
        let full = format_configuration(&config);
        assert_eq!(full, "1s2 2s2 2p6 3s2 3p6 4s2 3d6");
    }

    #[test]
    fn iron_short_config() {
        let config = electron_configuration(26).unwrap();
        let short = format_configuration_short(&config, 26);
        assert_eq!(short, "[Ar] 4s2 3d6");
    }

    #[test]
    fn chromium_exception() {
        // Cr (Z=24): [Ar] 3d5 4s1 (half-filled d shell stability)
        let config = electron_configuration(24).unwrap();
        let short = format_configuration_short(&config, 24);
        assert_eq!(short, "[Ar] 4s1 3d5");
    }

    #[test]
    fn copper_exception() {
        // Cu (Z=29): [Ar] 3d10 4s1 (filled d shell stability)
        let config = electron_configuration(29).unwrap();
        let short = format_configuration_short(&config, 29);
        assert_eq!(short, "[Ar] 4s1 3d10");
    }

    #[test]
    fn hydrogen_config() {
        let config = electron_configuration(1).unwrap();
        assert_eq!(config.len(), 1);
        assert_eq!(config[0].n, 1);
        assert_eq!(config[0].orbital, OrbitalType::S);
        assert_eq!(config[0].electrons, 1);
    }

    #[test]
    fn helium_config() {
        let config = electron_configuration(2).unwrap();
        let full = format_configuration(&config);
        assert_eq!(full, "1s2");
    }

    #[test]
    fn noble_gas_neon() {
        let config = electron_configuration(10).unwrap();
        let full = format_configuration(&config);
        assert_eq!(full, "1s2 2s2 2p6");
    }

    #[test]
    fn total_electrons_match_z() {
        for z in 1..=36 {
            let config = electron_configuration(z).unwrap();
            let total: u32 = config.iter().map(|e| e.electrons).sum();
            assert_eq!(total, z, "Z={z}: total electrons={total}");
        }
    }

    #[test]
    fn ionization_energy_hydrogen() {
        let ie = ionization_energy_ev(1).unwrap();
        assert!((ie - 13.598).abs() < 0.01);
    }

    #[test]
    fn ionization_energy_helium() {
        let ie = ionization_energy_ev(2).unwrap();
        assert!((ie - 24.587).abs() < 0.01);
    }

    #[test]
    fn ionization_energy_noble_gases_high() {
        // Noble gases should have the highest IE in their period
        let ne = ionization_energy_ev(10).unwrap();
        let na = ionization_energy_ev(11).unwrap();
        assert!(ne > na, "Ne IE={ne} should be > Na IE={na}");

        let ar = ionization_energy_ev(18).unwrap();
        let k = ionization_energy_ev(19).unwrap();
        assert!(ar > k, "Ar IE={ar} should be > K IE={k}");
    }

    #[test]
    fn ionization_energy_invalid() {
        assert!(ionization_energy_ev(0).is_err());
        assert!(ionization_energy_ev(37).is_err());
    }

    #[test]
    fn quantum_numbers_valid() {
        assert!(QuantumNumbers::new(1, 0, 0, 1).is_ok());
        assert!(QuantumNumbers::new(2, 1, -1, -1).is_ok());
        assert!(QuantumNumbers::new(3, 2, 2, 1).is_ok());
    }

    #[test]
    fn quantum_numbers_invalid() {
        assert!(QuantumNumbers::new(0, 0, 0, 1).is_err()); // n=0
        assert!(QuantumNumbers::new(1, 1, 0, 1).is_err()); // l >= n
        assert!(QuantumNumbers::new(2, 1, 2, 1).is_err()); // |ml| > l
        assert!(QuantumNumbers::new(1, 0, 0, 2).is_err()); // ms not +/-1
    }

    #[test]
    fn serde_roundtrip_orbital_type() {
        let o = OrbitalType::D;
        let json = serde_json::to_string(&o).unwrap();
        let back: OrbitalType = serde_json::from_str(&json).unwrap();
        assert_eq!(o, back);
    }

    #[test]
    fn serde_roundtrip_quantum_numbers() {
        let qn = QuantumNumbers::new(3, 2, -1, 1).unwrap();
        let json = serde_json::to_string(&qn).unwrap();
        let back: QuantumNumbers = serde_json::from_str(&json).unwrap();
        assert_eq!(qn, back);
    }

    #[test]
    fn serde_roundtrip_orbital_filling() {
        let of = OrbitalFilling {
            n: 3,
            orbital: OrbitalType::D,
            electrons: 6,
        };
        let json = serde_json::to_string(&of).unwrap();
        let back: OrbitalFilling = serde_json::from_str(&json).unwrap();
        assert_eq!(of, back);
    }

    #[test]
    fn orbital_max_electrons() {
        assert_eq!(OrbitalType::S.max_electrons(), 2);
        assert_eq!(OrbitalType::P.max_electrons(), 6);
        assert_eq!(OrbitalType::D.max_electrons(), 10);
        assert_eq!(OrbitalType::F.max_electrons(), 14);
    }
}
