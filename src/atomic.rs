//! Atomic structure: quantum numbers, electron configurations, and spectral lines.
//!
//! Implements the hydrogen-like spectral line formula (Rydberg), electron
//! configurations via the Aufbau principle with Madelung's rule, and
//! ionization energies from NIST data.

use crate::constants::{
    ATOMIC_UNIT_TIME_S, BOHR_MAGNETON_EV_T, BOHR_RADIUS, C, ELECTRON_ANOMALOUS_MOMENT,
    ELECTRON_MASS_MEV, FINE_STRUCTURE, H_EV_S, HC_EV_NM, PROTON_MASS_MEV, RYDBERG, RYDBERG_EV,
};
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
/// Follows the Aufbau principle with Madelung's rule, including all 20
/// ground-state exceptions listed by NIST (Cr, Cu, Nb, Mo, Ru, Rh, Pd, Ag, La,
/// Ce, Gd, Pt, Au, Ac, Th, Pa, U, Np, Cm, Lr).
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

    // Remove zero-electron entries (from exceptions that empty a subshell)
    config.retain(|e| e.electrons > 0);

    Ok(config)
}

/// Sets the electron count for a specific subshell in the configuration.
/// If the subshell doesn't exist in the config, it is inserted at the end.
fn set_subshell(config: &mut Vec<OrbitalFilling>, n: u32, orbital: OrbitalType, electrons: u32) {
    for entry in config.iter_mut() {
        if entry.n == n && entry.orbital == orbital {
            entry.electrons = electrons;
            return;
        }
    }
    // Subshell not found — add it (needed for La, Ce, Gd, Ac, Th, Pa, U, Np, Cm, Lr)
    if electrons > 0 {
        config.push(OrbitalFilling {
            n,
            orbital,
            electrons,
        });
    }
}

/// Applies known ground-state electron configuration exceptions (NIST).
///
/// Source: NIST Atomic Spectra Database ground-state configurations: the 20
/// Aufbau/Madelung exceptions for Z ≤ 103. For Z ≥ 104 relativistic
/// calculations predict Madelung-order ground states (e.g. Ds 6d⁸7s²,
/// Rg 6d⁹7s²; Smits et al., Phys. Rep. 1035, 1 (2023)), so no exceptions apply.
#[allow(clippy::too_many_lines)]
fn apply_exceptions(z: u32, config: &mut Vec<OrbitalFilling>) {
    match z {
        // --- Period 4 (3d block) ---
        24 => {
            // Cr: [Ar] 3d5 4s1 (half-filled d shell)
            set_subshell(config, 4, OrbitalType::S, 1);
            set_subshell(config, 3, OrbitalType::D, 5);
        }
        29 => {
            // Cu: [Ar] 3d10 4s1 (filled d shell)
            set_subshell(config, 4, OrbitalType::S, 1);
            set_subshell(config, 3, OrbitalType::D, 10);
        }

        // --- Period 5 (4d block) ---
        41 => {
            // Nb: [Kr] 4d4 5s1
            set_subshell(config, 5, OrbitalType::S, 1);
            set_subshell(config, 4, OrbitalType::D, 4);
        }
        42 => {
            // Mo: [Kr] 4d5 5s1 (half-filled d shell)
            set_subshell(config, 5, OrbitalType::S, 1);
            set_subshell(config, 4, OrbitalType::D, 5);
        }
        44 => {
            // Ru: [Kr] 4d7 5s1
            set_subshell(config, 5, OrbitalType::S, 1);
            set_subshell(config, 4, OrbitalType::D, 7);
        }
        45 => {
            // Rh: [Kr] 4d8 5s1
            set_subshell(config, 5, OrbitalType::S, 1);
            set_subshell(config, 4, OrbitalType::D, 8);
        }
        46 => {
            // Pd: [Kr] 4d10 5s0 (filled d, empty s)
            set_subshell(config, 5, OrbitalType::S, 0);
            set_subshell(config, 4, OrbitalType::D, 10);
        }
        47 => {
            // Ag: [Kr] 4d10 5s1 (filled d shell)
            set_subshell(config, 5, OrbitalType::S, 1);
            set_subshell(config, 4, OrbitalType::D, 10);
        }

        // --- Period 6: Lanthanides (4f block) ---
        57 => {
            // La: [Xe] 5d1 6s2 (electron goes to 5d, not 4f)
            set_subshell(config, 4, OrbitalType::F, 0);
            set_subshell(config, 5, OrbitalType::D, 1);
        }
        58 => {
            // Ce: [Xe] 4f1 5d1 6s2
            set_subshell(config, 4, OrbitalType::F, 1);
            set_subshell(config, 5, OrbitalType::D, 1);
        }
        64 => {
            // Gd: [Xe] 4f7 5d1 6s2 (half-filled f shell)
            set_subshell(config, 4, OrbitalType::F, 7);
            set_subshell(config, 5, OrbitalType::D, 1);
        }

        // --- Period 6: 5d block ---
        78 => {
            // Pt: [Xe] 4f14 5d9 6s1
            set_subshell(config, 6, OrbitalType::S, 1);
            set_subshell(config, 5, OrbitalType::D, 9);
        }
        79 => {
            // Au: [Xe] 4f14 5d10 6s1 (filled d shell)
            set_subshell(config, 6, OrbitalType::S, 1);
            set_subshell(config, 5, OrbitalType::D, 10);
        }

        // --- Period 7: Actinides (5f block) ---
        89 => {
            // Ac: [Rn] 6d1 7s2 (electron goes to 6d, not 5f)
            set_subshell(config, 5, OrbitalType::F, 0);
            set_subshell(config, 6, OrbitalType::D, 1);
        }
        90 => {
            // Th: [Rn] 6d2 7s2 (both electrons go to 6d)
            set_subshell(config, 5, OrbitalType::F, 0);
            set_subshell(config, 6, OrbitalType::D, 2);
        }
        91 => {
            // Pa: [Rn] 5f2 6d1 7s2
            set_subshell(config, 5, OrbitalType::F, 2);
            set_subshell(config, 6, OrbitalType::D, 1);
        }
        92 => {
            // U: [Rn] 5f3 6d1 7s2
            set_subshell(config, 5, OrbitalType::F, 3);
            set_subshell(config, 6, OrbitalType::D, 1);
        }
        93 => {
            // Np: [Rn] 5f4 6d1 7s2
            set_subshell(config, 5, OrbitalType::F, 4);
            set_subshell(config, 6, OrbitalType::D, 1);
        }
        96 => {
            // Cm: [Rn] 5f7 6d1 7s2 (half-filled f shell)
            set_subshell(config, 5, OrbitalType::F, 7);
            set_subshell(config, 6, OrbitalType::D, 1);
        }

        // --- Period 7: 6d block / superheavy ---
        103 => {
            // Lr: [Rn] 5f14 7s2 7p1 (relativistic: 7p lower than 6d)
            set_subshell(config, 6, OrbitalType::D, 0);
            set_subshell(config, 7, OrbitalType::P, 1);
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
/// Subshells are listed in filling order, e.g. `"[Ar] 4s2 3d6"` for iron (Z=26).
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

/// Calculates the vacuum wavelength of a spectral line in nanometers using
/// the Rydberg formula for a hydrogen-like atom with an **infinitely heavy**
/// nucleus.
///
/// 1/λ = R∞ Z² |1/n1² − 1/n2²|
///
/// The finite nuclear mass shifts real lines to longer wavelength by the
/// factor 1 + m_e/M (5.4e-4 for hydrogen: Hα 656.112 nm here, 656.470 nm
/// with the proton mass); use [`spectral_line_vacuum_nm`] for that. Air
/// wavelengths (Hα 656.28 nm) are shorter by the refractive index of air.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] if n1 or n2 is 0, or n1 == n2,
/// and [`TanmatraError::InvalidAtomicNumber`] if z is 0.
#[inline]
pub fn spectral_line_nm(z: u32, n1: u32, n2: u32) -> Result<f64, TanmatraError> {
    if z == 0 {
        return Err(TanmatraError::InvalidAtomicNumber(z));
    }
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

/// Returns the reduced-mass factor μ/m_e = 1/(1 + m_e/M) for a hydrogen-like
/// ion of nuclear charge `z` and mass number `a`.
///
/// The nuclear mass comes from the AME2020 atomic mass when tabulated
/// (minus Z electron masses), otherwise from the semi-empirical mass formula;
/// either is far more accurate than the 1e-7 relative precision this factor
/// needs.
///
/// # Errors
///
/// Returns an error if (z, a) is not a valid nucleus.
pub fn reduced_mass_factor(z: u32, a: u32) -> Result<f64, TanmatraError> {
    let nucleus = crate::nucleus::Nucleus::new(z, a)?;
    let nuclear_mass_mev = if z == 1 && a == 1 {
        PROTON_MASS_MEV
    } else {
        nucleus.experimental_atomic_mass_amu().map_or_else(
            || nucleus.nuclear_mass(),
            |m| m * crate::constants::AMU_MEV - z as f64 * ELECTRON_MASS_MEV,
        )
    };
    Ok(1.0 / (1.0 + ELECTRON_MASS_MEV / nuclear_mass_mev))
}

/// Vacuum wavelength of a hydrogen-like line including the nuclear
/// reduced-mass correction, in nanometers.
///
/// λ = λ∞ / (μ/m_e), with λ∞ from [`spectral_line_nm`].
/// Hydrogen (z=1, a=1) Hα: 656.470 nm.
///
/// # Errors
///
/// Returns an error for invalid quantum numbers or an invalid nucleus.
pub fn spectral_line_vacuum_nm(z: u32, a: u32, n1: u32, n2: u32) -> Result<f64, TanmatraError> {
    Ok(spectral_line_nm(z, n1, n2)? / reduced_mass_factor(z, a)?)
}

/// Calculates the energy of a hydrogen-like level with the first-order
/// fine-structure correction (infinite nuclear mass).
///
/// E_nj = −R∞hc Z² / n² × [1 + (αZ)²/n² × (n/(j+1/2) − 3/4)]
///
/// where R∞hc = 13.605693122990 eV (CODATA 2022), α is the fine-structure
/// constant and j = l ± 1/2 is passed as 2j. For the exact Dirac energy use
/// [`dirac_binding_energy_ev`].
///
/// Returns the energy in eV (negative, bound state).
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] if n = 0 or 2j is not an
/// odd number in 1..=2n−1.
pub fn hydrogen_level_energy_ev(z: u32, n: u32, two_j: u32) -> Result<f64, TanmatraError> {
    if n == 0 {
        return Err(TanmatraError::InvalidQuantumNumbers(String::from(
            "n must be >= 1",
        )));
    }
    if two_j.is_multiple_of(2) || two_j > 2 * n - 1 {
        return Err(TanmatraError::InvalidQuantumNumbers(alloc::format!(
            "2j={two_j} invalid for n={n}: must be odd and <= 2n-1"
        )));
    }

    let z_f = z as f64;
    let n_f = n as f64;
    let j_plus_half = f64::midpoint(two_j as f64, 1.0);

    let e0 = -RYDBERG_EV * z_f * z_f / (n_f * n_f);
    let az = FINE_STRUCTURE * z_f;
    let correction = 1.0 + az * az / (n_f * n_f) * (n_f / j_plus_half - 0.75);

    Ok(e0 * correction)
}

/// Calculates the wavelength of a spectral line with fine-structure correction.
///
/// Uses the first-order fine-structure energies of
/// [`hydrogen_level_energy_ev`] (infinite nuclear mass, vacuum wavelength).
///
/// Parameters:
/// - `z`: atomic number
/// - `n1`, `two_j1`: lower level (principal quantum number, 2*j)
/// - `n2`, `two_j2`: upper level (principal quantum number, 2*j)
///
/// Returns wavelength in nanometers.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] if any quantum numbers are invalid.
pub fn spectral_line_fine_nm(
    z: u32,
    n1: u32,
    two_j1: u32,
    n2: u32,
    two_j2: u32,
) -> Result<f64, TanmatraError> {
    let e1 = hydrogen_level_energy_ev(z, n1, two_j1)?;
    let e2 = hydrogen_level_energy_ev(z, n2, two_j2)?;

    let delta_e = (e2 - e1).abs();
    if delta_e < 1e-30 {
        return Err(TanmatraError::InvalidQuantumNumbers(String::from(
            "transition energy is zero",
        )));
    }

    // λ = hc/ΔE
    Ok(HC_EV_NM / delta_e)
}

// ---------------------------------------------------------------------------
// Named spectral series
// ---------------------------------------------------------------------------

/// Returns the Lyman series wavelengths (n_upper -> n=1) in nanometers.
///
/// The Lyman series lies in the ultraviolet (91.2 - 121.6 nm).
/// Returns lines for n_upper = 2 through `n_max`.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] if `n_max` < 2.
pub fn lyman_series(z: u32, n_max: u32) -> Result<Vec<(u32, f64)>, TanmatraError> {
    if n_max < 2 {
        return Err(TanmatraError::InvalidQuantumNumbers(String::from(
            "Lyman series requires n_max >= 2",
        )));
    }
    let mut lines = Vec::new();
    for n in 2..=n_max {
        lines.push((n, spectral_line_nm(z, 1, n)?));
    }
    Ok(lines)
}

/// Returns the Balmer series wavelengths (n_upper -> n=2) in nanometers.
///
/// The Balmer series spans the visible and near-UV (364.6 - 656.3 nm).
/// Returns lines for n_upper = 3 through `n_max`.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] if `n_max` < 3.
pub fn balmer_series(z: u32, n_max: u32) -> Result<Vec<(u32, f64)>, TanmatraError> {
    if n_max < 3 {
        return Err(TanmatraError::InvalidQuantumNumbers(String::from(
            "Balmer series requires n_max >= 3",
        )));
    }
    let mut lines = Vec::new();
    for n in 3..=n_max {
        lines.push((n, spectral_line_nm(z, 2, n)?));
    }
    Ok(lines)
}

/// Returns the Paschen series wavelengths (n_upper -> n=3) in nanometers.
///
/// The Paschen series lies in the near-infrared (820.4 - 1875 nm).
/// Returns lines for n_upper = 4 through `n_max`.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] if `n_max` < 4.
pub fn paschen_series(z: u32, n_max: u32) -> Result<Vec<(u32, f64)>, TanmatraError> {
    if n_max < 4 {
        return Err(TanmatraError::InvalidQuantumNumbers(String::from(
            "Paschen series requires n_max >= 4",
        )));
    }
    let mut lines = Vec::new();
    for n in 4..=n_max {
        lines.push((n, spectral_line_nm(z, 3, n)?));
    }
    Ok(lines)
}

/// Returns the Brackett series wavelengths (n_upper -> n=4) in nanometers.
///
/// The Brackett series lies in the infrared (1458 - 4051 nm).
/// Returns lines for n_upper = 5 through `n_max`.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] if `n_max` < 5.
pub fn brackett_series(z: u32, n_max: u32) -> Result<Vec<(u32, f64)>, TanmatraError> {
    if n_max < 5 {
        return Err(TanmatraError::InvalidQuantumNumbers(String::from(
            "Brackett series requires n_max >= 5",
        )));
    }
    let mut lines = Vec::new();
    for n in 5..=n_max {
        lines.push((n, spectral_line_nm(z, 4, n)?));
    }
    Ok(lines)
}

/// Returns the Pfund series wavelengths (n_upper -> n=5) in nanometers.
///
/// The Pfund series lies in the far-infrared (2279 - 7460 nm).
/// Returns lines for n_upper = 6 through `n_max`.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] if `n_max` < 6.
pub fn pfund_series(z: u32, n_max: u32) -> Result<Vec<(u32, f64)>, TanmatraError> {
    if n_max < 6 {
        return Err(TanmatraError::InvalidQuantumNumbers(String::from(
            "Pfund series requires n_max >= 6",
        )));
    }
    let mut lines = Vec::new();
    for n in 6..=n_max {
        lines.push((n, spectral_line_nm(z, 5, n)?));
    }
    Ok(lines)
}

// ---------------------------------------------------------------------------
// Zeeman and Stark effects
// ---------------------------------------------------------------------------

/// Calculates the Lande g-factor for a single-electron level.
///
/// g_J = 1 + [J(J+1) + S(S+1) − L(L+1)] / [2J(J+1)]
///
/// with S = 1/2, L = l, J = j and g_s = 2.
///
/// Parameters: `l` (orbital), `two_j` (2*total angular momentum).
/// Returns 0.0 if j is not l ± 1/2 (no such level).
#[must_use]
#[inline]
pub fn lande_g_factor(l: u32, two_j: u32) -> f64 {
    if two_j + 1 != 2 * l + 2 && two_j + 1 != 2 * l {
        return 0.0;
    }
    let j = two_j as f64 / 2.0;
    let l_f = l as f64;
    let s = 0.5;

    let j_j1 = j * (j + 1.0);
    if j_j1 < 1e-30 {
        return 0.0;
    }

    1.0 + (j_j1 + s * (s + 1.0) - l_f * (l_f + 1.0)) / (2.0 * j_j1)
}

/// Calculates the anomalous Zeeman energy splitting in eV.
///
/// ΔE = m_j * g_J * μ_B * B
///
/// where m_j ranges from -j to +j in integer steps.
///
/// Parameters:
/// - `l`: orbital quantum number
/// - `two_j`: 2*j (total angular momentum)
/// - `two_mj`: 2*m_j (magnetic quantum number projection)
/// - `b_tesla`: magnetic field strength in Tesla
///
/// Returns the energy shift in eV.
#[must_use]
#[inline]
pub fn zeeman_splitting_ev(l: u32, two_j: u32, two_mj: i32, b_tesla: f64) -> f64 {
    let mj = two_mj as f64 / 2.0;
    let g = lande_g_factor(l, two_j);
    mj * g * BOHR_MAGNETON_EV_T * b_tesla
}

/// Calculates the linear Stark effect energy shift for hydrogen (Z = 1).
///
/// ΔE = (3/2) n k e a₀ F
///
/// where k = n₁ − n₂ is the difference of parabolic quantum numbers
/// (n₁ + n₂ + |m| + 1 = n, so |k| ≤ n − 1) and F the field strength.
///
/// Parameters:
/// - `n`: principal quantum number
/// - `parabolic_index`: k = n₁ − n₂, from −(n−1) to (n−1)
/// - `e_field_v_per_m`: electric field strength in V/m
///
/// Returns the energy shift in eV, or 0.0 if |k| > n − 1.
#[must_use]
#[inline]
pub fn stark_shift_hydrogen_ev(n: u32, parabolic_index: i32, e_field_v_per_m: f64) -> f64 {
    if n == 0 || parabolic_index.unsigned_abs() > n - 1 {
        return 0.0;
    }
    // e·a₀·F in joules divided by e gives eV: a₀ [m] × F [V/m].
    1.5 * n as f64 * parabolic_index as f64 * BOHR_RADIUS * e_field_v_per_m
}

// ---------------------------------------------------------------------------
// QED corrections
// ---------------------------------------------------------------------------

/// Bethe logarithms ln k₀(n, l) for hydrogen (Drake & Swainson,
/// Phys. Rev. A 41, 1243 (1990)).
fn bethe_log(n: u32, l: u32) -> Option<f64> {
    Some(match (n, l) {
        (1, 0) => 2.984_128_556,
        (2, 0) => 2.811_769_893,
        (2, 1) => -0.030_016_709,
        (3, 0) => 2.767_663_612,
        (3, 1) => -0.038_190_229,
        (3, 2) => -0.005_232_148,
        (4, 0) => 2.749_811_840,
        (4, 1) => -0.041_954_895,
        (4, 2) => -0.006_740_939,
        (4, 3) => -0.001_733_661,
        _ => return None,
    })
}

/// One-loop QED (Lamb) shift of the hydrogen-like level (n, l, j) in eV,
/// infinite nuclear mass.
///
/// ΔE = (α/π) (Zα)⁴ m_e c² / n³ × [A₄₁ ln(Zα)⁻² + A₄₀ + Zα A₅₀]
///
/// Self-energy: A₄₁ = 4/3 δ_l0, A₄₀ = −(4/3) ln k₀(n,l) + 10/9 δ_l0
/// − (1 − δ_l0)/(2κ(2l+1)), A₅₀ = (139/32 − 2 ln 2) π δ_l0.
/// Vacuum polarization (Uehling): A₄₀ = −4/15 δ_l0, A₅₀ = (5/48) π δ_l0.
/// κ = −(l+1) for j = l + 1/2 and κ = l for j = l − 1/2.
/// (Eides, Grotch & Shelyuto, Phys. Rep. 342, 63 (2001); Mohr et al.,
/// CODATA 2018 Rev. Mod. Phys. 93, 025010 (2021).)
///
/// The Zα expansion is accurate to ≈0.2% for hydrogen (2S½−2P½:
/// 1059.6 MHz vs 1057.845 MHz measured) and ≈3% for He⁺; it is not valid for
/// high Z. Bethe logarithms are tabulated for n ≤ 4.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] for n > 4, l ≥ n, or j ≠ l ± 1/2.
pub fn lamb_shift_nlj_ev(z: u32, n: u32, l: u32, two_j: u32) -> Result<f64, TanmatraError> {
    if n == 0 || l >= n || (two_j != 2 * l + 1 && two_j + 1 != 2 * l) {
        return Err(TanmatraError::InvalidQuantumNumbers(alloc::format!(
            "invalid level n={n} l={l} 2j={two_j}"
        )));
    }
    let ln_k0 = bethe_log(n, l).ok_or_else(|| {
        TanmatraError::InvalidQuantumNumbers(alloc::format!(
            "Bethe logarithm not tabulated for n={n} l={l}"
        ))
    })?;
    let za = z as f64 * FINE_STRUCTURE;
    let nf = n as f64;
    let prefactor = FINE_STRUCTURE / core::f64::consts::PI * za.powi(4) * ELECTRON_MASS_MEV * 1e6
        / (nf * nf * nf);
    let pi = core::f64::consts::PI;
    let coefficient = if l == 0 {
        let a41 = 4.0 / 3.0 * libm::log(1.0 / (za * za));
        let a40 = -4.0 / 3.0 * ln_k0 + 10.0 / 9.0 - 4.0 / 15.0;
        let a50 = (139.0 / 32.0 - 2.0 * core::f64::consts::LN_2) * pi + 5.0 / 48.0 * pi;
        a41 + a40 + za * a50
    } else {
        let kappa = if two_j == 2 * l + 1 {
            -(l as f64 + 1.0)
        } else {
            l as f64
        };
        -4.0 / 3.0 * ln_k0 - 1.0 / (2.0 * kappa * (2.0 * l as f64 + 1.0))
    };
    Ok(prefactor * coefficient)
}

/// One-loop QED (Lamb) shift of a hydrogen-like (n, l) level in eV.
///
/// For s states this is the shift of nS½ from [`lamb_shift_nlj_ev`]. For
/// l > 0 it is the (2j+1)-weighted mean over j = l ± 1/2, for which the
/// κ-dependent terms cancel. Returns 0.0 where [`lamb_shift_nlj_ev`] would
/// return an error (n = 0, l ≥ n, or n > 4).
#[must_use]
pub fn lamb_shift_ev(z: u32, n: u32, l: u32) -> f64 {
    if l == 0 {
        return lamb_shift_nlj_ev(z, n, 0, 1).unwrap_or(0.0);
    }
    let lo = lamb_shift_nlj_ev(z, n, l, 2 * l - 1).unwrap_or(0.0);
    let hi = lamb_shift_nlj_ev(z, n, l, 2 * l + 1).unwrap_or(0.0);
    (2.0 * l as f64 * lo + (2.0 * l as f64 + 2.0) * hi) / (4.0 * l as f64 + 2.0)
}

/// Vacuum polarization (Uehling potential) shift of a hydrogen-like s level
/// in eV, infinite nuclear mass.
///
/// ΔE_VP = (α/π)(Zα)⁴ m_e c² / n³ × [−4/15 + (5π/48) Zα]
///
/// For hydrogen 2S: −26.89 MHz (leading term −27.13 MHz, α(Zα)⁵ term
/// +0.24 MHz). Returns 0.0 for l > 0, where these terms vanish.
#[must_use]
#[inline]
pub fn vacuum_polarization_ev(z: u32, n: u32, l: u32) -> f64 {
    if n == 0 || l >= n || l != 0 {
        return 0.0;
    }
    let za = z as f64 * FINE_STRUCTURE;
    let nf = n as f64;
    FINE_STRUCTURE / core::f64::consts::PI * za.powi(4) * ELECTRON_MASS_MEV * 1e6 / (nf * nf * nf)
        * (-4.0 / 15.0 + 5.0 * core::f64::consts::PI / 48.0 * za)
}

// ---------------------------------------------------------------------------
// Hydrogen wavefunctions
// ---------------------------------------------------------------------------

/// Coefficients c_i of the polynomial P(ρ) = ρ^l L_{n−l−1}^{(2l+1)}(ρ) in powers
/// of ρ (index = power), and the normalization N of
/// R_nl = N P(ρ) e^{−ρ/2}, ρ = 2Zr/(n a₀), in units of a₀^{−3/2}.
fn hydrogen_radial_polynomial(z: u32, n: u32, l: u32) -> (Vec<f64>, f64) {
    let k = n - l - 1;
    let alpha = 2 * l + 1;
    let mut coeffs = alloc::vec![0.0; (l + k + 1) as usize];
    // L_k^α(x) = Σ_i (−1)^i C(k+α, k−i) x^i / i!
    for i in 0..=k {
        let mut binom = 1.0_f64; // C(k+α, k−i)
        let top = k + alpha;
        let choose = k - i;
        for t in 0..choose {
            binom *= f64::from(top - t) / f64::from(t + 1);
        }
        let mut fact = 1.0_f64;
        for t in 1..=i {
            fact *= f64::from(t);
        }
        let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
        coeffs[(l + i) as usize] = sign * binom / fact;
    }
    // N = sqrt((2Z/n)³ (n−l−1)! / (2n (n+l)!))
    let zf = f64::from(z);
    let nf = f64::from(n);
    let mut ratio = 1.0_f64; // (n−l−1)!/(n+l)! = 1/Π_{i=n−l}^{n+l} i
    for t in (n - l)..=(n + l) {
        ratio /= f64::from(t);
    }
    let scale = 2.0 * zf / nf;
    let norm = libm::sqrt(scale * scale * scale * ratio / (2.0 * nf));
    (coeffs, norm)
}

/// Computes the radial wavefunction R_nl(r) for a hydrogen-like atom
/// (infinite nuclear mass).
///
/// R_nl(r) = N ρ^l e^{−ρ/2} L_{n−l−1}^{(2l+1)}(ρ),  ρ = 2Zr/(n a₀),
/// N = √((2Z/n a₀)³ (n−l−1)! / (2n (n+l)!)),
///
/// with the generalized Laguerre polynomial L, normalized so that
/// ∫₀^∞ R² r² dr = 1.
///
/// Parameters:
/// - `z`: atomic number (≥ 1)
/// - `n`: principal quantum number (1..=60)
/// - `l`: orbital angular momentum (0 to n−1)
/// - `r_bohr`: radial distance in units of the Bohr radius (r/a₀)
///
/// Returns R_nl(r) in units of a₀^(−3/2).
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] if n = 0, n > 60, l ≥ n or z = 0.
pub fn radial_wavefunction(z: u32, n: u32, l: u32, r_bohr: f64) -> Result<f64, TanmatraError> {
    if n == 0 || n > 60 {
        return Err(TanmatraError::InvalidQuantumNumbers(alloc::format!(
            "n={n} not supported (1-60)"
        )));
    }
    if l >= n {
        return Err(TanmatraError::InvalidQuantumNumbers(alloc::format!(
            "l={l} must be < n={n}"
        )));
    }
    if z == 0 {
        return Err(TanmatraError::InvalidAtomicNumber(z));
    }
    let (coeffs, norm) = hydrogen_radial_polynomial(z, n, l);
    let rho = 2.0 * f64::from(z) * r_bohr / f64::from(n);
    let mut poly = 0.0;
    for c in coeffs.iter().rev() {
        poly = poly * rho + c;
    }
    Ok(norm * poly * libm::exp(-rho / 2.0))
}

/// Exact dipole radial integral ⟨n l | r | n′ l′⟩ = ∫ R_nl R_n′l′ r³ dr in units
/// of a₀/Z (hydrogen-like, infinite nuclear mass).
fn radial_dipole_integral(n: u32, l: u32, n2: u32, l2: u32) -> f64 {
    let (c1, norm1) = hydrogen_radial_polynomial(1, n, l);
    let (c2, norm2) = hydrogen_radial_polynomial(1, n2, l2);
    // R1 R2 r³ = N1 N2 Σ_i Σ_j c1_i c2_j (2/n)^i (2/n2)^j r^(i+j+3) e^{−(1/n + 1/n2) r}
    let s1 = 2.0 / f64::from(n);
    let s2 = 2.0 / f64::from(n2);
    let decay = 1.0 / f64::from(n) + 1.0 / f64::from(n2);
    let mut total = 0.0;
    for (i, a) in c1.iter().enumerate() {
        if *a == 0.0 {
            continue;
        }
        for (jdx, b) in c2.iter().enumerate() {
            if *b == 0.0 {
                continue;
            }
            let power = i + jdx + 3;
            // ∫ r^p e^{−c r} dr = p! / c^{p+1}, accumulated as a product to avoid overflow.
            let mut integral = 1.0 / decay;
            for t in 1..=power {
                integral *= t as f64 / decay;
            }
            total += a * b * libm::pow(s1, i as f64) * libm::pow(s2, jdx as f64) * integral;
        }
    }
    norm1 * norm2 * total
}

/// Computes the radial probability density |R_nl(r)|² * r² for hydrogen-like atoms.
///
/// This is the probability of finding the electron at distance r (per unit r).
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] if quantum numbers are invalid.
#[inline]
pub fn radial_probability_density(
    z: u32,
    n: u32,
    l: u32,
    r_bohr: f64,
) -> Result<f64, TanmatraError> {
    let rnl = radial_wavefunction(z, n, l, r_bohr)?;
    Ok(rnl * rnl * r_bohr * r_bohr)
}

// ---------------------------------------------------------------------------
// Selection rules and transition probabilities
// ---------------------------------------------------------------------------

/// Result of checking electric dipole selection rules for a transition.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum TransitionType {
    /// Electric dipole (E1) — allowed transition.
    ElectricDipole,
    /// Forbidden by electric dipole selection rules.
    Forbidden,
}

/// Checks if a transition satisfies electric dipole (E1) selection rules.
///
/// For hydrogen-like atoms, E1 selection rules require:
/// - Δl = ±1
/// - Δm_l = 0, ±1
/// - Δm_s = 0
/// - Parity change (automatic from Δl = ±1)
///
/// Parameters: `l1`, `l2` (orbital quantum numbers of initial/final states).
#[must_use]
#[inline]
pub fn check_selection_rules(l1: u32, l2: u32) -> TransitionType {
    let delta_l = (l1 as i64 - l2 as i64).unsigned_abs();
    if delta_l == 1 {
        TransitionType::ElectricDipole
    } else {
        TransitionType::Forbidden
    }
}

/// Checks if a transition satisfies the full E1 selection rules including m_l.
///
/// Rules: Δl = ±1, Δm_l = 0, ±1.
#[must_use]
#[inline]
pub fn check_selection_rules_full(l1: u32, ml1: i32, l2: u32, ml2: i32) -> TransitionType {
    let delta_l = (l1 as i64 - l2 as i64).unsigned_abs();
    let delta_ml = (ml1 - ml2).unsigned_abs();
    if delta_l == 1 && delta_ml <= 1 {
        TransitionType::ElectricDipole
    } else {
        TransitionType::Forbidden
    }
}

/// Calculates the Einstein A coefficient (spontaneous emission rate) for an
/// electric dipole transition (n_u, l_u) → (n_l, l_l) in a hydrogen-like ion,
/// summed over final and averaged over initial magnetic substates (infinite
/// nuclear mass, non-relativistic).
///
/// A = (4/3) α³ ω³ × max(l_u, l_l)/(2l_u + 1) × |⟨n_l l_l|r|n_u l_u⟩|²  (atomic units)
///
/// with ω = Z²(1/(2n_l²) − 1/(2n_u²)) E_h and the radial integral evaluated
/// exactly from the hydrogenic wavefunctions (∝ 1/Z), so A ∝ Z⁴. Converted to
/// s⁻¹ with the atomic unit of time ħ/E_h. For a finite nuclear mass multiply
/// by μ/m_e (see [`reduced_mass_factor`]): H Lyα gives 6.2684e8 s⁻¹ here and
/// 6.2649e8 s⁻¹ with μ, matching NIST ASD.
///
/// Returns the rate in s⁻¹.
///
/// # Errors
///
/// Returns an error if the transition violates Δl = ±1, if n_u ≤ n_l, if
/// l ≥ n for either level, if n_u > 60, or if z = 0.
pub fn einstein_a_coefficient(
    z: u32,
    n_upper: u32,
    l_upper: u32,
    n_lower: u32,
    l_lower: u32,
) -> Result<f64, TanmatraError> {
    if check_selection_rules(l_upper, l_lower) == TransitionType::Forbidden {
        return Err(TanmatraError::InvalidQuantumNumbers(alloc::format!(
            "transition ({n_upper},{l_upper})->({n_lower},{l_lower}) is E1-forbidden"
        )));
    }
    if n_upper <= n_lower {
        return Err(TanmatraError::InvalidQuantumNumbers(String::from(
            "upper level must have larger n",
        )));
    }
    if n_lower == 0 || l_upper >= n_upper || l_lower >= n_lower || n_upper > 60 {
        return Err(TanmatraError::InvalidQuantumNumbers(alloc::format!(
            "invalid levels ({n_upper},{l_upper})->({n_lower},{l_lower})"
        )));
    }
    if z == 0 {
        return Err(TanmatraError::InvalidAtomicNumber(z));
    }

    let zf = f64::from(z);
    let nu = f64::from(n_upper);
    let nl = f64::from(n_lower);
    let omega_au = zf * zf * 0.5 * (1.0 / (nl * nl) - 1.0 / (nu * nu));
    let radial_au = radial_dipole_integral(n_upper, l_upper, n_lower, l_lower) / zf;
    let l_max = f64::from(l_upper.max(l_lower));
    let angular = l_max / f64::from(2 * l_upper + 1);
    let alpha3 = FINE_STRUCTURE * FINE_STRUCTURE * FINE_STRUCTURE;

    let rate_au =
        4.0 / 3.0 * alpha3 * omega_au * omega_au * omega_au * angular * radial_au * radial_au;
    Ok(rate_au / ATOMIC_UNIT_TIME_S)
}

/// Calculates the Einstein B coefficient for stimulated emission.
///
/// B₂₁ = A₂₁ c³ / (8π h ν³)
///
/// where ν is the transition frequency (infinite nuclear mass). This B refers
/// to the spectral energy density per unit frequency; absorption follows from
/// g₁B₁₂ = g₂B₂₁.
///
/// Returns B in m³/(J·s²).
///
/// # Errors
///
/// Returns error if the transition violates selection rules.
pub fn einstein_b_coefficient(
    z: u32,
    n_upper: u32,
    l_upper: u32,
    n_lower: u32,
    l_lower: u32,
) -> Result<f64, TanmatraError> {
    let a21 = einstein_a_coefficient(z, n_upper, l_upper, n_lower, l_lower)?;

    let zf = z as f64;
    let n1 = n_lower as f64;
    let n2 = n_upper as f64;

    let delta_e_ev = RYDBERG_EV * zf * zf * (1.0 / (n1 * n1) - 1.0 / (n2 * n2));
    let freq_hz = delta_e_ev / H_EV_S;
    let h_joule_s = H_EV_S * crate::constants::ELEMENTARY_CHARGE;

    if freq_hz <= 0.0 {
        return Ok(0.0);
    }

    Ok(a21 * C * C * C / (8.0 * core::f64::consts::PI * h_joule_s * freq_hz.powi(3)))
}

// ---------------------------------------------------------------------------
// Electron affinities
// ---------------------------------------------------------------------------

/// Electron affinity of a neutral atom.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum ElectronAffinity {
    /// The anion is bound; the electron affinity in eV.
    Bound(f64),
    /// No bound anion exists (negative electron affinity).
    Unbound,
    /// No measurement or accepted calculation is available.
    Unknown,
}

/// Electron affinities for Z=1..=118.
///
/// Measured values: Andersen, Haugen & Hotop, J. Phys. Chem. Ref. Data 28,
/// 1511 (1999), superseded where newer laser-photodetachment measurements exist
/// (cited per entry). Elements without measurements use the cited calculations
/// or estimates; superheavy values are from Smits et al., Phys. Rep. 1035, 1
/// (2023), Table 2, unless noted.
#[allow(clippy::too_many_lines)]
const ELECTRON_AFFINITY_EV: [ElectronAffinity; 118] = [
    ElectronAffinity::Bound(0.754195),    // H (Z=1) Lykke 1991; AHH99
    ElectronAffinity::Unbound,            // He (Z=2) AHH99 <0
    ElectronAffinity::Bound(0.618049),    // Li (Z=3) AHH99
    ElectronAffinity::Unbound,            // Be (Z=4) AHH99 <0
    ElectronAffinity::Bound(0.279723),    // B (Z=5) AHH99
    ElectronAffinity::Bound(1.2621226),   // C (Z=6) Bresteau 2016 (AHH99 1.262118(20))
    ElectronAffinity::Unbound,            // N (Z=7) AHH99 -0.07(2)
    ElectronAffinity::Bound(1.461112972), // O (Z=8) Kristiansson 2022
    ElectronAffinity::Bound(3.4011887),   // F (Z=9) AHH99
    ElectronAffinity::Unbound,            // Ne (Z=10) AHH99 <0
    ElectronAffinity::Bound(0.547926),    // Na (Z=11) AHH99
    ElectronAffinity::Unbound,            // Mg (Z=12) AHH99 <0
    ElectronAffinity::Bound(0.43283),     // Al (Z=13) AHH99
    ElectronAffinity::Bound(1.389521),    // Si (Z=14) AHH99
    ElectronAffinity::Bound(0.746609),    // P (Z=15) Pelaez 2011 (AHH99 0.7465(3))
    ElectronAffinity::Bound(2.0771029),   // S (Z=16) AHH99
    ElectronAffinity::Bound(3.612724),    // Cl (Z=17) AHH99
    ElectronAffinity::Unbound,            // Ar (Z=18) AHH99 <0
    ElectronAffinity::Bound(0.501459),    // K (Z=19) AHH99
    ElectronAffinity::Bound(0.02455),     // Ca (Z=20) AHH99
    ElectronAffinity::Bound(0.179378),    // Sc (Z=21) Lu et al. JCP 2023 (AHH99 0.188(20))
    ElectronAffinity::Bound(0.07554),     // Ti (Z=22) Tang 2018 (AHH99 0.084(9))
    ElectronAffinity::Bound(0.52766),     // V (Z=23) Fu 2016 (AHH99 0.525(12))
    ElectronAffinity::Bound(0.67584),     // Cr (Z=24) Bilodeau 1998 5451.0(10) cm-1; AHH99
    ElectronAffinity::Unbound,            // Mn (Z=25) AHH99 <0
    ElectronAffinity::Bound(0.153236),    // Fe (Z=26) Chen 2016 (AHH99 0.151(3))
    ElectronAffinity::Bound(0.662256),    // Co (Z=27) Chen&Ning 2016 (AHH99 0.6633(6))
    ElectronAffinity::Bound(1.15716),     // Ni (Z=28) Scheer 1998 1157.16(12) meV; AHH99
    ElectronAffinity::Bound(1.23578),     // Cu (Z=29) Bilodeau 1998; AHH99
    ElectronAffinity::Unbound,            // Zn (Z=30) AHH99 <0
    ElectronAffinity::Bound(0.301166),    // Ga (Z=31) Tang 2020 (AHH99 0.41(4))
    ElectronAffinity::Bound(1.2326764),   // Ge (Z=32) Bresteau 2015 (AHH99 1.232712(15))
    ElectronAffinity::Bound(0.804486),    // As (Z=33) Blondel&Drag 2025 (AHH99 0.814(8))
    ElectronAffinity::Bound(2.020667),    // Se (Z=34) Zhang 2026 (AHH99 2.02067(2))
    ElectronAffinity::Bound(3.363588),    // Br (Z=35) AHH99
    ElectronAffinity::Unbound,            // Kr (Z=36) AHH99 <0
    ElectronAffinity::Bound(0.485916),    // Rb (Z=37) AHH99
    ElectronAffinity::Bound(0.05206),     // Sr (Z=38) Andersen 1997 52.06(6) meV; AHH99
    ElectronAffinity::Bound(0.31129),     // Y (Z=39) Lu et al. JCP 2023 (AHH99 0.307(12))
    ElectronAffinity::Bound(0.433283),    // Zr (Z=40) Fu 2017 (AHH99 0.426(14))
    ElectronAffinity::Bound(0.9174),      // Nb (Z=41) Luo 2016 (AHH99 0.893(25))
    ElectronAffinity::Bound(0.7472),      // Mo (Z=42) Bilodeau 1998 6027(2) cm-1; AHH99
    ElectronAffinity::Bound(0.55),        // Tc (Z=43) AHH99 semi-empirical est.
    ElectronAffinity::Bound(1.04638),     // Ru (Z=44) AHH99 (NingLu2022 1.04627(2), W)
    ElectronAffinity::Bound(1.14289),     // Rh (Z=45) Scheer 1998; AHH99
    ElectronAffinity::Bound(0.56214),     // Pd (Z=46) Scheer 1998; AHH99
    ElectronAffinity::Bound(1.30447),     // Ag (Z=47) AHH99 10521.3(2) cm-1 (Bilodeau 1998)
    ElectronAffinity::Unbound,            // Cd (Z=48) AHH99 <0
    ElectronAffinity::Bound(0.38392),     // In (Z=49) Walter 2010 (AHH99 0.404(9))
    ElectronAffinity::Bound(1.11207),     // Sn (Z=50) Vandevraye 2013 (AHH99 1.112066(15))
    ElectronAffinity::Bound(1.047401),    // Sb (Z=51) Scheer 1997 8447.86(15) cm-1; AHH99
    ElectronAffinity::Bound(1.970875),    // Te (Z=52) AHH99
    ElectronAffinity::Bound(3.0590465),   // I (Z=53) Pelaez 2009 (AHH99 3.059038(10))
    ElectronAffinity::Unbound,            // Xe (Z=54) AHH99 <0
    ElectronAffinity::Bound(0.4715983),   // Cs (Z=55) Navarro Navarrete 2024 (AHH99 0.471626(25))
    ElectronAffinity::Bound(0.14462),     // Ba (Z=56) AHH99
    ElectronAffinity::Bound(0.557546),    // La (Z=57) Blondel 2020 / Lu 2019 (AHH99 0.47(2))
    ElectronAffinity::Bound(0.60016),     // Ce (Z=58) Fu 2020 (AHH99: <0.5 est.)
    ElectronAffinity::Bound(0.10923),     // Pr (Z=59) Fu 2020 PRA 101 022502
    ElectronAffinity::Bound(0.09748),     // Nd (Z=60) Fu 2020 PRA 101 022502
    ElectronAffinity::Bound(0.129),       // Pm (Z=61) Felfli 2009 calc (uncertain)
    ElectronAffinity::Bound(0.162),       // Sm (Z=62) Felfli 2009 calc (uncertain)
    ElectronAffinity::Bound(0.116),       // Eu (Z=63) Cheng&Castleman 2015
    ElectronAffinity::Bound(0.212),       // Gd (Z=64) NingLu 2022 (primary not verified)
    ElectronAffinity::Bound(0.13131),     // Tb (Z=65) Fu 2020 PRA 101 022502
    ElectronAffinity::Bound(0.015),       // Dy (Z=66) Nadeau 1997 AMS (uncertain)
    ElectronAffinity::Bound(0.338),       // Ho (Z=67) Felfli 2009 calc (uncertain)
    ElectronAffinity::Bound(0.312),       // Er (Z=68) Felfli 2009 calc (uncertain)
    ElectronAffinity::Bound(1.029),       // Tm (Z=69) Davis&Thompson 2001
    ElectronAffinity::Unbound,            // Yb (Z=70) CRC est. -0.02
    ElectronAffinity::Bound(0.23882),     // Lu (Z=71) Fu 2019
    ElectronAffinity::Bound(0.178),       // Hf (Z=72) Tang 2018 (AHH99 ~0)
    ElectronAffinity::Bound(0.322),       // Ta (Z=73) AHH99 (NingLu2022 0.328859(23), W)
    ElectronAffinity::Bound(0.815),       // W (Z=74) AHH99 (NingLu2022 0.816500(82), W)
    ElectronAffinity::Bound(0.060396),    // Re (Z=75) Chen&Ning 2017 (AHH99 0.15(15))
    ElectronAffinity::Bound(1.0778),      // Os (Z=76) AHH99 (NingLu2022 1.077661(24), W)
    ElectronAffinity::Bound(1.564057),    // Ir (Z=77) Lu 2020 (Bilodeau 1999/AHH99 1.56436(15))
    ElectronAffinity::Bound(2.1251),      // Pt (Z=78) Bilodeau 1999 17140.1(4) cm-1; AHH99
    ElectronAffinity::Bound(2.30861),     // Au (Z=79) AHH99
    ElectronAffinity::Unbound,            // Hg (Z=80) AHH99 <0
    ElectronAffinity::Bound(0.320053),    // Tl (Z=81) Walter 2020 (AHH99 0.377(13))
    ElectronAffinity::Bound(0.356721),    // Pb (Z=82) Bresteau 2019 (AHH99 0.364(8))
    ElectronAffinity::Bound(0.942363),    // Bi (Z=83) AHH99
    ElectronAffinity::Bound(1.4), // Po (Z=84) Li 2012 MCDHF calc; AHH99 SE est 1.9(3) (unmeasured)
    ElectronAffinity::Bound(2.41578), // At (Z=85) Leimbach 2020
    ElectronAffinity::Unbound,    // Rn (Z=86) AHH99 <0
    ElectronAffinity::Bound(0.491), // Fr (Z=87) Landau 2001 calc
    ElectronAffinity::Bound(0.1), // Ra (Z=88) CRC/Andersen est. (uncertain)
    ElectronAffinity::Bound(0.35), // Ac (Z=89) CRC est. (uncertain)
    ElectronAffinity::Bound(0.60769), // Th (Z=90) Tang 2019 PRL
    ElectronAffinity::Bound(0.55), // Pa (Z=91) est. (uncertain)
    ElectronAffinity::Bound(0.31497), // U (Z=92) Tang 2021 PRA; Ciborowski 2021 0.309(25)
    ElectronAffinity::Bound(0.48), // Np (Z=93) est. (uncertain)
    ElectronAffinity::Unbound,    // Pu (Z=94) est. -0.50 (uncertain)
    ElectronAffinity::Bound(0.1), // Am (Z=95) est. (uncertain)
    ElectronAffinity::Bound(0.28), // Cm (Z=96) est. (uncertain)
    ElectronAffinity::Unbound,    // Bk (Z=97) est. -1.72 (uncertain)
    ElectronAffinity::Unbound,    // Cf (Z=98) est. -1.01 (uncertain)
    ElectronAffinity::Unbound,    // Es (Z=99) est. -0.30 (uncertain)
    ElectronAffinity::Bound(0.35), // Fm (Z=100) est. (uncertain)
    ElectronAffinity::Bound(0.98), // Md (Z=101) est. (uncertain)
    ElectronAffinity::Unbound,    // No (Z=102) est. -2.33 (uncertain)
    ElectronAffinity::Bound(0.446), // Lr (Z=103) Guo 2024 PRA 110 022817 (FSCC 2007: 0.476)
    ElectronAffinity::Unknown,    // Rf (Z=104) no value in Smits 2023 Tab.2
    ElectronAffinity::Bound(1.189), // Db (Z=105) Smits 2023 Tab.2
    ElectronAffinity::Unknown,    // Sg (Z=106) no value in Smits 2023
    ElectronAffinity::Unknown,    // Bh (Z=107) no value
    ElectronAffinity::Unknown,    // Hs (Z=108) no value
    ElectronAffinity::Unknown,    // Mt (Z=109) no value
    ElectronAffinity::Bound(0.83), // Ds (Z=110) Smits 2023 Tab.2
    ElectronAffinity::Bound(1.56), // Rg (Z=111) Eliav 1994 abstract 1.56 (Eliav 2015 1.565; Smits 2023 1.97)
    ElectronAffinity::Unbound,     // Cn (Z=112) Borschevsky slides / Smits 2023: 0
    ElectronAffinity::Bound(0.776), // Nh (Z=113) Guo 2022 JPB (Borschevsky 0.69; Smits 0.73)
    ElectronAffinity::Unbound,     // Fl (Z=114) Borschevsky slides / Smits 2023: 0 (no EA)
    ElectronAffinity::Bound(0.313), // Mc (Z=115) Smits 2023 (Borschevsky 0.366)
    ElectronAffinity::Bound(0.776), // Lv (Z=116) Smits 2023
    ElectronAffinity::Bound(1.602), // Ts (Z=117) Smits 2023 (Borschevsky 1.719)
    ElectronAffinity::Bound(0.076), // Og (Z=118) Kaygorodov 2021 (Eliav 1996: 0.056)
];

/// Returns the electron affinity of the neutral atom, distinguishing bound,
/// unbound and unknown anions.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidAtomicNumber`] if Z is 0 or > 118.
#[inline]
pub fn electron_affinity(z: u32) -> Result<ElectronAffinity, TanmatraError> {
    if z == 0 || z > 118 {
        return Err(TanmatraError::InvalidAtomicNumber(z));
    }
    Ok(ELECTRON_AFFINITY_EV[(z - 1) as usize])
}

/// Returns the electron affinity in eV for the given atomic number.
///
/// The electron affinity is the energy released when an electron is added
/// to a neutral atom: X + e⁻ → X⁻ + EA.
///
/// Returns 0.0 both when no bound anion exists and when no value is known;
/// use [`electron_affinity`] to tell these apart.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidAtomicNumber`] if Z is 0 or > 118.
#[inline]
pub fn electron_affinity_ev(z: u32) -> Result<f64, TanmatraError> {
    Ok(match electron_affinity(z)? {
        ElectronAffinity::Bound(ev) => ev,
        ElectronAffinity::Unbound | ElectronAffinity::Unknown => 0.0,
    })
}

/// First ionization energies for Z=1 to Z=118 in eV.
///
/// Z=1–108: NIST Atomic Spectra Database, Ionization Energies Data
/// (Kramida, Ralchenko, Reader & NIST ASD Team), retrieved 2026-09. Values NIST
/// marks as theoretical, interpolated or estimated are noted per entry.
/// <https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html>
///
/// Z=109–118: relativistic calculations compiled in Smits et al.,
/// Phys. Rep. 1035, 1 (2023), Table 2 (no measurements exist).
#[allow(clippy::too_many_lines)]
const IONIZATION_ENERGIES_EV: [f64; 118] = [
    13.598434599702, // H (Z=1)
    24.587389011,    // He (Z=2)
    5.391714996,     // Li (Z=3)
    9.322699,        // Be (Z=4)
    8.298019,        // B (Z=5)
    11.260288,       // C (Z=6)
    14.53413,        // N (Z=7)
    13.618055,       // O (Z=8)
    17.42282,        // F (Z=9)
    21.564541,       // Ne (Z=10)
    5.13907696,      // Na (Z=11)
    7.646236,        // Mg (Z=12)
    5.985769,        // Al (Z=13)
    8.15168,         // Si (Z=14)
    10.486686,       // P (Z=15)
    10.3600167,      // S (Z=16)
    12.967633,       // Cl (Z=17)
    15.7596119,      // Ar (Z=18)
    4.34066373,      // K (Z=19)
    6.113154921,     // Ca (Z=20)
    6.56149,         // Sc (Z=21)
    6.82812,         // Ti (Z=22)
    6.746187,        // V (Z=23)
    6.76651,         // Cr (Z=24)
    7.434038,        // Mn (Z=25)
    7.9024681,       // Fe (Z=26)
    7.88101,         // Co (Z=27)
    7.639878,        // Ni (Z=28)
    7.72638,         // Cu (Z=29)
    9.394197,        // Zn (Z=30)
    5.999302,        // Ga (Z=31)
    7.899435,        // Ge (Z=32)
    9.78855,         // As (Z=33)
    9.752368,        // Se (Z=34)
    11.81381,        // Br (Z=35)
    13.9996055,      // Kr (Z=36)
    4.1771281,       // Rb (Z=37)
    5.69486745,      // Sr (Z=38)
    6.21726,         // Y (Z=39)
    6.634126,        // Zr (Z=40)
    6.75885,         // Nb (Z=41)
    7.09243,         // Mo (Z=42)
    7.11938,         // Tc (Z=43)
    7.3605,          // Ru (Z=44)
    7.4589,          // Rh (Z=45)
    8.336839,        // Pd (Z=46)
    7.576234,        // Ag (Z=47)
    8.99382,         // Cd (Z=48)
    5.7863558,       // In (Z=49)
    7.343918,        // Sn (Z=50)
    8.608389,        // Sb (Z=51)
    9.009808,        // Te (Z=52)
    10.451236,       // I (Z=53)
    12.1298437,      // Xe (Z=54)
    3.89390572743,   // Cs (Z=55)
    5.2116646,       // Ba (Z=56)
    5.5769,          // La (Z=57)
    5.5386,          // Ce (Z=58)
    5.4702,          // Pr (Z=59) — NIST ASD
    5.52475,         // Nd (Z=60)
    5.58187,         // Pm (Z=61)
    5.643722,        // Sm (Z=62)
    5.670385,        // Eu (Z=63)
    6.1498,          // Gd (Z=64)
    5.8638,          // Tb (Z=65)
    5.939061,        // Dy (Z=66)
    6.0215,          // Ho (Z=67)
    6.1077,          // Er (Z=68)
    6.184402,        // Tm (Z=69)
    6.25416,         // Yb (Z=70)
    5.425871,        // Lu (Z=71)
    6.82507,         // Hf (Z=72)
    7.549571,        // Ta (Z=73)
    7.86403,         // W (Z=74)
    7.83352,         // Re (Z=75)
    8.43823,         // Os (Z=76)
    8.96702,         // Ir (Z=77)
    8.95883,         // Pt (Z=78)
    9.225554,        // Au (Z=79)
    10.437504,       // Hg (Z=80)
    6.1082873,       // Tl (Z=81)
    7.4166799,       // Pb (Z=82)
    7.285516,        // Bi (Z=83)
    8.41807,         // Po (Z=84)
    9.31751,         // At (Z=85)
    10.7485,         // Rn (Z=86)
    4.0727411,       // Fr (Z=87)
    5.2784239,       // Ra (Z=88)
    5.380235,        // Ac (Z=89)
    6.3067,          // Th (Z=90)
    5.89,            // Pa (Z=91) — NIST ASD (theory/estimate)
    6.19405,         // U (Z=92)
    6.265608,        // Np (Z=93)
    6.02576,         // Pu (Z=94)
    5.97381,         // Am (Z=95)
    5.992241,        // Cm (Z=96)
    6.19785,         // Bk (Z=97)
    6.281878,        // Cf (Z=98)
    6.3684,          // Es (Z=99)
    6.5,             // Fm (Z=100) — NIST ASD (theory/estimate)
    6.58,            // Md (Z=101) — NIST ASD (theory/estimate)
    6.62621,         // No (Z=102)
    4.96,            // Lr (Z=103)
    6.02,            // Rf (Z=104) — NIST ASD (theory/estimate)
    6.8,             // Db (Z=105) — NIST ASD (theory/estimate)
    7.8,             // Sg (Z=106) — NIST ASD (theory/estimate)
    7.7,             // Bh (Z=107) — NIST ASD (theory/estimate)
    7.6,             // Hs (Z=108) — NIST ASD (theory/estimate)
    10.4,            // Mt (Z=109) — Smits et al. 2023, theory
    9.562,           // Ds (Z=110) — Smits et al. 2023, theory
    11.03,           // Rg (Z=111) — Smits et al. 2023, theory
    12.02,           // Cn (Z=112) — Smits et al. 2023, theory
    7.49,            // Nh (Z=113) — Smits et al. 2023, theory
    8.65,            // Fl (Z=114) — Smits et al. 2023, theory
    5.574,           // Mc (Z=115) — Smits et al. 2023, theory
    6.855,           // Lv (Z=116) — Smits et al. 2023, theory
    7.654,           // Ts (Z=117) — Smits et al. 2023, theory
    8.888,           // Og (Z=118) — Smits et al. 2023, theory
];

/// Returns the first ionization energy in eV for the given atomic number.
///
/// Covers Z=1 (hydrogen) through Z=118 (oganesson).
/// Z=1-103 from NIST ASD experimental values.
/// Z=104-118 from relativistic theoretical predictions.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidAtomicNumber`] if Z is 0 or > 118.
#[inline]
pub fn ionization_energy_ev(z: u32) -> Result<f64, TanmatraError> {
    if z == 0 || z > 118 {
        return Err(TanmatraError::InvalidAtomicNumber(z));
    }
    Ok(IONIZATION_ENERGIES_EV[(z - 1) as usize])
}

// ---------------------------------------------------------------------------
// Relativistic quantum: Dirac equation, hyperfine, anomalous moment, Breit
// ---------------------------------------------------------------------------

/// Calculates the exact Dirac energy for a hydrogen-like atom in MeV.
///
/// Uses the exact solution of the Dirac equation:
///
/// E_nj = m_e c² / √(1 + [αZ / (n - δ)]²)
///
/// where δ = j + 1/2 - √((j+1/2)² - (αZ)²).
///
/// Parameters:
/// - `z`: atomic number
/// - `n`: principal quantum number (>= 1)
/// - `two_j`: twice the total angular momentum quantum number (must be odd, >= 1)
///
/// Returns the total Dirac energy in MeV (rest mass + binding).
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] if quantum numbers are invalid.
pub fn dirac_energy_mev(z: u32, n: u32, two_j: u32) -> Result<f64, TanmatraError> {
    if n == 0 {
        return Err(TanmatraError::InvalidQuantumNumbers(String::from(
            "n must be >= 1",
        )));
    }
    if two_j == 0 || two_j.is_multiple_of(2) || two_j > 2 * n - 1 {
        return Err(TanmatraError::InvalidQuantumNumbers(alloc::format!(
            "two_j={two_j} must be odd, >= 1 and <= 2n-1 (half-integer j <= n - 1/2)"
        )));
    }

    let alpha = FINE_STRUCTURE;
    let me_c2 = ELECTRON_MASS_MEV;
    let zf = z as f64;
    let nf = n as f64;
    let j_plus_half = f64::midpoint(two_j as f64, 1.0); // j + 1/2

    let az = alpha * zf;
    let az2 = az * az;

    // δ = j + 1/2 - √((j+1/2)² - (αZ)²)
    let discriminant = j_plus_half * j_plus_half - az2;
    if discriminant <= 0.0 {
        return Err(TanmatraError::InvalidQuantumNumbers(alloc::format!(
            "αZ={az} too large for j+1/2={j_plus_half}: supercritical"
        )));
    }
    let delta = j_plus_half - libm::sqrt(discriminant);

    // n_eff = n - δ
    let n_eff = nf - delta;
    if n_eff <= 0.0 {
        return Err(TanmatraError::InvalidQuantumNumbers(alloc::format!(
            "effective quantum number n-δ={n_eff} must be positive for n={n}, two_j={two_j}"
        )));
    }

    // E = m_e c² / √(1 + (αZ / n_eff)²)
    let ratio = az / n_eff;
    let energy = me_c2 / libm::sqrt(1.0 + ratio * ratio);

    Ok(energy)
}

/// Calculates the Dirac binding energy for a hydrogen-like atom in eV.
///
/// Returns m_e c² - E_nj (positive for bound states), where E_nj is the
/// exact Dirac energy.
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] if quantum numbers are invalid.
#[inline]
pub fn dirac_binding_energy_ev(z: u32, n: u32, two_j: u32) -> Result<f64, TanmatraError> {
    let e_dirac = dirac_energy_mev(z, n, two_j)?;
    let me_c2 = ELECTRON_MASS_MEV;
    Ok((me_c2 - e_dirac) * 1e6)
}

/// Calculates the total relativistic correction to hydrogen-like energy levels in eV.
///
/// Returns the difference between the exact Dirac binding energy and the
/// non-relativistic binding energy: ΔE = E_Dirac_binding - E_nonrel_binding.
///
/// This captures all orders of relativistic corrections (fine structure,
/// second-order, third-order, etc.).
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] if quantum numbers are invalid.
pub fn relativistic_correction_ev(z: u32, n: u32, two_j: u32) -> Result<f64, TanmatraError> {
    let binding_dirac = dirac_binding_energy_ev(z, n, two_j)?;

    // Non-relativistic binding energy ½α²m_ec² Z²/n², from the same constants
    // as the Dirac energy so the difference is self-consistent.
    let zf = z as f64;
    let nf = n as f64;
    let rydberg_ev = 0.5 * FINE_STRUCTURE * FINE_STRUCTURE * ELECTRON_MASS_MEV * 1e6;
    let binding_nonrel = rydberg_ev * zf * zf / (nf * nf);

    Ok(binding_dirac - binding_nonrel)
}

/// Calculates the magnetic-dipole hyperfine splitting of an s state of a
/// hydrogen-like ion with a nuclear spin-1/2 nucleus, in eV.
///
/// ΔE = (4/3) α⁴ Z³ g_I (m_e/m_p) m_e c² / n³  (Fermi contact interaction,
/// leading order: no reduced-mass, QED or relativistic corrections)
///
/// For hydrogen 1s this gives 5.8776e-6 eV (1421.2 MHz; measured 1420.406 MHz).
/// For nuclear spin I ≠ 1/2 use [`hyperfine_splitting_spin_ev`].
///
/// Parameters:
/// - `z`: atomic number
/// - `n`: principal quantum number
/// - `nuclear_g_factor`: g_I = μ_I/(I μ_N) (5.5856946893 for the proton)
#[must_use]
#[inline]
pub fn hyperfine_splitting_ev(z: u32, n: u32, nuclear_g_factor: f64) -> f64 {
    hyperfine_splitting_spin_ev(z, n, nuclear_g_factor, 0.5)
}

/// Calculates the magnetic-dipole hyperfine splitting between F = I + 1/2 and
/// F = I − 1/2 of an s state of a hydrogen-like ion, in eV.
///
/// ΔE = A (I + 1/2),  A = (4/3) α⁴ Z³ g_I (m_e/m_p) m_e c² / n³
///
/// Leading order (Fermi contact term). Deuterium 1s (I = 1,
/// g_d = 0.8574382335): 327.23 MHz (measured 327.384 MHz).
///
/// Returns 0.0 for n = 0 or I ≤ 0.
#[must_use]
#[inline]
pub fn hyperfine_splitting_spin_ev(
    z: u32,
    n: u32,
    nuclear_g_factor: f64,
    nuclear_spin: f64,
) -> f64 {
    if n == 0 || nuclear_spin <= 0.0 {
        return 0.0;
    }
    let zf = z as f64;
    let nf = n as f64;
    let alpha4 = FINE_STRUCTURE * FINE_STRUCTURE * FINE_STRUCTURE * FINE_STRUCTURE;
    let mass_ratio = ELECTRON_MASS_MEV / PROTON_MASS_MEV;
    let a_const = (4.0 / 3.0)
        * alpha4
        * zf
        * zf
        * zf
        * mass_ratio
        * nuclear_g_factor
        * ELECTRON_MASS_MEV
        * 1e6
        / (nf * nf * nf);
    a_const * (nuclear_spin + 0.5)
}

/// Returns the free-electron g-factor including QED corrections.
///
/// g_e = 2(1 + a_e) where a_e is the anomalous magnetic moment.
///
/// CODATA 2022: |g_e| = 2.00231930436092.
#[must_use]
#[inline]
pub fn electron_g_factor() -> f64 {
    2.0 * (1.0 + ELECTRON_ANOMALOUS_MOMENT)
}

/// Calculates the bound-electron g-factor of the 1s ground state of a
/// hydrogen-like ion (Breit 1928 Dirac value times the free-electron QED factor).
///
/// g_bound = (2/3)(1 + 2√(1 - (αZ)²)) × (1 + a_e)
///
/// This gives the leading relativistic correction to the electron g-factor
/// in the Coulomb field of the nucleus.
///
/// For Z=1: g_bound ≈ 2.002 (very close to free g_e).
/// For high Z: g_bound decreases due to relativistic effects.
#[must_use]
#[inline]
pub fn bound_electron_g_factor(z: u32) -> f64 {
    let alpha = FINE_STRUCTURE;
    let az = alpha * z as f64;
    let az2 = az * az;

    if az2 >= 1.0 {
        // Supercritical: formula invalid, return 0 rather than panic
        return 0.0;
    }

    let breit = (2.0 / 3.0) * (1.0 + 2.0 * libm::sqrt(1.0 - az2));
    breit * (1.0 + ELECTRON_ANOMALOUS_MOMENT)
}

/// Calculates the strong-field (Paschen–Back) Zeeman shift of a
/// single-electron state in eV, using the QED electron g-factor.
///
/// ΔE = (m_l + g_e m_s) μ_B B,  g_e = 2.00231930436092 (CODATA 2022)
///
/// Valid when the Zeeman energy greatly exceeds the fine-structure splitting,
/// so that l and s decouple. For weak fields use [`zeeman_splitting_ev`]
/// (m_j g_J μ_B B).
///
/// Parameters:
/// - `b_tesla`: magnetic field strength in Tesla
/// - `ml`: orbital magnetic quantum number
/// - `ms`: spin magnetic quantum number (+1 or -1, representing ±1/2)
///
/// Returns the energy shift in eV.
#[must_use]
#[inline]
pub fn anomalous_zeeman_splitting_ev(b_tesla: f64, ml: i32, ms: i32) -> f64 {
    let g_e = electron_g_factor();
    let ms_half = ms as f64 / 2.0;
    (ml as f64 + g_e * ms_half) * BOHR_MAGNETON_EV_T * b_tesla
}

/// Calculates the leading-order Breit (magnetic) interaction energy of the
/// 1s² ¹S₀ ground state of a helium-like ion in eV.
///
/// In the Breit–Pauli reduction the two-electron magnetic terms for 1s²
/// ¹S₀ reduce to the spin–spin contact term, since orbit–orbit and
/// spin–other-orbit terms vanish for two s electrons in a spin singlet
/// (Bethe & Salpeter, Quantum Mechanics of One- and Two-Electron Atoms, §38–39):
///
/// ΔE = −(8π/3) α² ⟨s₁·s₂⟩ ⟨δ³(r₁₂)⟩ = 2π α² ⟨δ³(r₁₂)⟩ E_h a₀³,
///
/// and with unscreened hydrogenic 1s orbitals ⟨δ³(r₁₂)⟩ = Z³/(8π a₀³):
///
/// ΔE = α² Z³ / 4 E_h  (positive: the singlet is raised).
///
/// This is the leading term of the 1/Z expansion. Electron correlation reduces
/// it substantially at low Z (for He the correlated value of the contact term
/// is ≈ 0.67 α² E_h, versus 2 α² E_h here) and relativistic corrections of
/// order (Zα)² increase it at high Z; the formula is most useful for
/// intermediate Z (roughly 10–40).
#[must_use]
#[inline]
pub fn breit_interaction_ev(z: u32) -> f64 {
    let zf = z as f64;
    let hartree_ev = 2.0 * RYDBERG_EV;
    FINE_STRUCTURE * FINE_STRUCTURE * zf * zf * zf / 4.0 * hartree_ev
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
        for z in 1..=118 {
            let config = electron_configuration(z).unwrap();
            let total: u32 = config.iter().map(|e| e.electrons).sum();
            assert_eq!(total, z, "Z={z}: total electrons={total}");
        }
    }

    #[test]
    fn niobium_exception() {
        // Nb (Z=41): [Kr] 4d4 5s1
        let config = electron_configuration(41).unwrap();
        let short = format_configuration_short(&config, 41);
        assert_eq!(short, "[Kr] 5s1 4d4");
    }

    #[test]
    fn palladium_exception() {
        // Pd (Z=46): [Kr] 4d10 (no 5s electrons at all)
        let config = electron_configuration(46).unwrap();
        let short = format_configuration_short(&config, 46);
        assert_eq!(short, "[Kr] 4d10");
    }

    #[test]
    fn lanthanum_exception() {
        // La (Z=57): [Xe] 5d1 6s2
        let config = electron_configuration(57).unwrap();
        let short = format_configuration_short(&config, 57);
        assert_eq!(short, "[Xe] 6s2 5d1");
    }

    #[test]
    fn gadolinium_exception() {
        // Gd (Z=64): [Xe] 4f7 5d1 6s2
        let config = electron_configuration(64).unwrap();
        let short = format_configuration_short(&config, 64);
        assert_eq!(short, "[Xe] 6s2 4f7 5d1");
    }

    #[test]
    fn thorium_exception() {
        // Th (Z=90): [Rn] 6d2 7s2
        let config = electron_configuration(90).unwrap();
        let short = format_configuration_short(&config, 90);
        assert_eq!(short, "[Rn] 7s2 6d2");
    }

    #[test]
    fn uranium_exception() {
        // U (Z=92): [Rn] 5f3 6d1 7s2
        let config = electron_configuration(92).unwrap();
        let short = format_configuration_short(&config, 92);
        assert_eq!(short, "[Rn] 7s2 5f3 6d1");
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
        assert!(ionization_energy_ev(119).is_err());
    }

    #[test]
    fn ionization_energy_cesium_lowest_alkali() {
        // Cs has the lowest IE of any non-superheavy element
        let cs = ionization_energy_ev(55).unwrap();
        assert!(cs < 4.0, "Cs IE={cs} should be < 4 eV");
    }

    #[test]
    fn ionization_energy_noble_gas_trend() {
        // Noble gases: He, Ne, Ar, Kr, Xe, Rn — should all be local maxima
        let rn = ionization_energy_ev(86).unwrap();
        let fr = ionization_energy_ev(87).unwrap();
        assert!(rn > fr, "Rn IE={rn} should be > Fr IE={fr}");
    }

    #[test]
    fn ionization_energy_lanthanide_range() {
        // Lanthanides (Z=57-71) should all be in ~5.4-6.3 eV range
        for z in 57..=71 {
            let ie = ionization_energy_ev(z).unwrap();
            assert!(ie > 5.0 && ie < 6.5, "Z={z} IE={ie} out of range");
        }
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

    // --- Fine-structure / relativistic correction tests ---

    #[test]
    fn hydrogen_ground_state_energy() {
        // H 1s1/2: should be approximately -13.6 eV
        let e = hydrogen_level_energy_ev(1, 1, 1).unwrap();
        assert!((e - (-13.6)).abs() < 0.01, "H 1s1/2 E={e} eV");
    }

    #[test]
    fn fine_structure_splits_levels() {
        // H n=2: 2s1/2 and 2p3/2 should have different energies
        // 2s1/2: 2j=1, 2p1/2: 2j=1, 2p3/2: 2j=3
        let e_s = hydrogen_level_energy_ev(1, 2, 1).unwrap(); // 2s1/2
        let e_p32 = hydrogen_level_energy_ev(1, 2, 3).unwrap(); // 2p3/2
        // Fine structure: 2p3/2 should have slightly different energy than 2s1/2
        assert!(
            (e_s - e_p32).abs() > 1e-6,
            "Fine structure should split n=2 levels"
        );
    }

    #[test]
    fn h_alpha_fine_structure() {
        // H-alpha with fine structure: 2p3/2 -> 3d5/2 transition
        // Should be close to the non-relativistic 656.3 nm
        let lambda = spectral_line_fine_nm(1, 2, 3, 3, 5).unwrap();
        assert!(
            (lambda - 656.3).abs() < 2.0,
            "H-alpha fine={lambda} nm, expected ~656 nm"
        );
    }

    #[test]
    fn fine_structure_invalid() {
        assert!(hydrogen_level_energy_ev(1, 0, 1).is_err());
        assert!(hydrogen_level_energy_ev(1, 1, 3).is_err());
    }

    // --- Zeeman/Stark tests ---

    #[test]
    fn lande_g_s_electron() {
        // 1s1/2 (l=0, j=1/2): g = 2.0
        let g = lande_g_factor(0, 1);
        assert!((g - 2.0).abs() < 1e-10, "g(1s1/2)={g}, expected 2.0");
    }

    #[test]
    fn lande_g_p32() {
        // 2p3/2 (l=1, j=3/2): g = 4/3
        let g = lande_g_factor(1, 3);
        assert!((g - 4.0 / 3.0).abs() < 1e-10, "g(2p3/2)={g}");
    }

    #[test]
    fn zeeman_zero_field() {
        let de = zeeman_splitting_ev(0, 1, 1, 0.0);
        assert!((de).abs() < 1e-30, "No splitting in zero field");
    }

    #[test]
    fn zeeman_proportional_to_field() {
        let de1 = zeeman_splitting_ev(0, 1, 1, 1.0);
        let de2 = zeeman_splitting_ev(0, 1, 1, 2.0);
        assert!(
            (de2 - 2.0 * de1).abs() < 1e-15,
            "Zeeman should be linear in B"
        );
    }

    #[test]
    fn stark_zero_field() {
        let de = stark_shift_hydrogen_ev(2, 1, 0.0);
        assert!((de).abs() < 1e-30);
    }

    #[test]
    fn stark_proportional_to_field() {
        let de1 = stark_shift_hydrogen_ev(2, 1, 1e6);
        let de2 = stark_shift_hydrogen_ev(2, 1, 2e6);
        assert!((de2 - 2.0 * de1).abs() < 1e-20);
    }

    #[test]
    fn stark_increases_with_n() {
        let de2 = stark_shift_hydrogen_ev(2, 1, 1e8).abs();
        let de5 = stark_shift_hydrogen_ev(5, 4, 1e8).abs();
        assert!(de5 > de2, "n=5 Stark should be larger than n=2");
    }

    // --- Electron affinity tests ---

    #[test]
    fn electron_affinity_fluorine_highest() {
        // F has the highest EA among period 2 elements
        let ea_f = electron_affinity_ev(9).unwrap();
        assert!(ea_f > 3.0, "F EA={ea_f}");
    }

    #[test]
    fn electron_affinity_chlorine() {
        let ea_cl = electron_affinity_ev(17).unwrap();
        assert!((ea_cl - 3.613).abs() < 0.01, "Cl EA={ea_cl}");
    }

    #[test]
    fn electron_affinity_noble_gases_zero() {
        for z in [2, 10, 18, 36, 54, 86] {
            let ea = electron_affinity_ev(z).unwrap();
            assert!((ea).abs() < 1e-10, "Noble gas Z={z} should have EA=0");
        }
    }

    #[test]
    fn electron_affinity_gold_high() {
        // Au has the highest EA among metals
        let ea_au = electron_affinity_ev(79).unwrap();
        assert!(ea_au > 2.0, "Au EA={ea_au}");
    }

    #[test]
    fn electron_affinity_invalid() {
        assert!(electron_affinity_ev(0).is_err());
        assert!(electron_affinity_ev(119).is_err());
    }

    // --- Wavefunction tests ---

    #[test]
    fn radial_1s_at_origin() {
        // R_10(0) = 2 * Z^{3/2} for hydrogen
        let r = radial_wavefunction(1, 1, 0, 0.0).unwrap();
        assert!((r - 2.0).abs() < 1e-10, "R_10(0)={r}, expected 2.0");
    }

    #[test]
    fn radial_1s_decays() {
        // R_10 should decrease with distance
        let r1 = radial_wavefunction(1, 1, 0, 1.0).unwrap().abs();
        let r5 = radial_wavefunction(1, 1, 0, 5.0).unwrap().abs();
        assert!(r1 > r5, "R_10 should decay with distance");
    }

    #[test]
    fn radial_2s_has_node() {
        // R_20 has a node at r = 2a0/Z = 2 for hydrogen
        let r = radial_wavefunction(1, 2, 0, 2.0).unwrap();
        assert!(r.abs() < 0.1, "R_20 should be near zero at r=2a0");
    }

    #[test]
    fn probability_density_positive() {
        let pd = radial_probability_density(1, 1, 0, 1.0).unwrap();
        assert!(pd > 0.0, "Probability density must be positive");
    }

    #[test]
    fn radial_invalid_quantum_numbers() {
        assert!(radial_wavefunction(1, 0, 0, 1.0).is_err()); // n=0
        assert!(radial_wavefunction(1, 1, 1, 1.0).is_err()); // l >= n
        assert!(radial_wavefunction(1, 61, 0, 1.0).is_err()); // n > 60
        assert!(radial_wavefunction(0, 1, 0, 1.0).is_err()); // z = 0
    }

    // --- Selection rules and Einstein coefficient tests ---

    #[test]
    fn selection_rules_allowed() {
        assert_eq!(check_selection_rules(0, 1), TransitionType::ElectricDipole);
        assert_eq!(check_selection_rules(1, 2), TransitionType::ElectricDipole);
        assert_eq!(check_selection_rules(2, 1), TransitionType::ElectricDipole);
    }

    #[test]
    fn selection_rules_forbidden() {
        assert_eq!(check_selection_rules(0, 0), TransitionType::Forbidden);
        assert_eq!(check_selection_rules(0, 2), TransitionType::Forbidden);
        assert_eq!(check_selection_rules(1, 3), TransitionType::Forbidden);
    }

    #[test]
    fn selection_rules_full_with_ml() {
        // Allowed: Δl=1, Δm_l=0
        assert_eq!(
            check_selection_rules_full(0, 0, 1, 0),
            TransitionType::ElectricDipole
        );
        // Forbidden: Δl=0
        assert_eq!(
            check_selection_rules_full(1, 0, 1, 0),
            TransitionType::Forbidden
        );
        // Forbidden: Δm_l=2
        assert_eq!(
            check_selection_rules_full(0, 0, 1, 2),
            TransitionType::Forbidden
        );
    }

    #[test]
    fn einstein_a_lyman_alpha() {
        // Lyman-alpha: 2p -> 1s, A ≈ 6.27e8 s⁻¹
        let a21 = einstein_a_coefficient(1, 2, 1, 1, 0).unwrap();
        assert!(a21 > 1e7, "Lyman-alpha A={a21}, should be ~6e8");
        assert!(a21 < 1e10, "Lyman-alpha A={a21} too large");
    }

    #[test]
    fn einstein_a_forbidden_transition() {
        // 2s -> 1s is forbidden (Δl = 0)
        assert!(einstein_a_coefficient(1, 2, 0, 1, 0).is_err());
    }

    #[test]
    fn einstein_b_positive() {
        let b21 = einstein_b_coefficient(1, 2, 1, 1, 0).unwrap();
        assert!(b21 > 0.0, "Einstein B should be positive");
    }

    #[test]
    fn serde_roundtrip_transition_type() {
        let tt = TransitionType::ElectricDipole;
        let json = serde_json::to_string(&tt).unwrap();
        let back: TransitionType = serde_json::from_str(&json).unwrap();
        assert_eq!(tt, back);
    }

    // --- QED correction tests ---

    #[test]
    fn lamb_shift_hydrogen_2s() {
        // Known: H 2S Lamb shift ≈ 4.37e-6 eV
        let shift = lamb_shift_ev(1, 2, 0);
        assert!((shift - 4.37e-6).abs() < 1e-7, "H 2S Lamb shift={shift} eV");
    }

    #[test]
    fn lamb_shift_grows_with_z() {
        assert!(lamb_shift_ev(2, 2, 0) > lamb_shift_ev(1, 2, 0));
    }

    #[test]
    fn lamb_shift_p_smaller_than_s() {
        let s_shift = lamb_shift_ev(1, 2, 0);
        let p_shift = lamb_shift_ev(1, 2, 1);
        assert!(
            s_shift > p_shift,
            "s-state Lamb shift should be larger than p-state"
        );
    }

    #[test]
    fn vacuum_polarization_negative() {
        let vp = vacuum_polarization_ev(1, 1, 0);
        assert!(vp < 0.0, "VP should be negative (downward shift)");
    }

    #[test]
    fn vacuum_polarization_small_vs_lamb() {
        let lamb = lamb_shift_ev(1, 2, 0);
        let vp = vacuum_polarization_ev(1, 2, 0).abs();
        assert!(vp < lamb, "VP should be smaller than Lamb shift");
    }

    // --- Named spectral series tests ---

    #[test]
    fn lyman_series_hydrogen() {
        let lines = lyman_series(1, 6).unwrap();
        assert_eq!(lines.len(), 5); // n=2..=6
        // Lyman-alpha: ~121.6 nm
        assert!((lines[0].1 - 121.6).abs() < 0.5, "Ly-α={}", lines[0].1);
        // Series limit: should converge toward 91.2 nm
        assert!(lines.last().unwrap().1 > 91.0);
    }

    #[test]
    fn balmer_series_hydrogen() {
        let lines = balmer_series(1, 8).unwrap();
        assert_eq!(lines.len(), 6); // n=3..=8
        // H-alpha: ~656.3 nm
        assert!((lines[0].1 - 656.3).abs() < 1.0, "Hα={}", lines[0].1);
        // H-beta: ~486.1 nm
        assert!((lines[1].1 - 486.1).abs() < 1.0, "Hβ={}", lines[1].1);
    }

    #[test]
    fn paschen_series_hydrogen() {
        let lines = paschen_series(1, 7).unwrap();
        assert_eq!(lines.len(), 4); // n=4..=7
        // Paschen-alpha: ~1875 nm
        assert!((lines[0].1 - 1875.0).abs() < 10.0, "Pa-α={}", lines[0].1);
    }

    #[test]
    fn brackett_series_hydrogen() {
        let lines = brackett_series(1, 8).unwrap();
        assert_eq!(lines.len(), 4); // n=5..=8
        // Brackett-alpha: ~4051 nm
        assert!((lines[0].1 - 4051.0).abs() < 20.0, "Br-α={}", lines[0].1);
    }

    #[test]
    fn pfund_series_hydrogen() {
        let lines = pfund_series(1, 9).unwrap();
        assert_eq!(lines.len(), 4); // n=6..=9
        // Pfund-alpha: ~7460 nm
        assert!((lines[0].1 - 7460.0).abs() < 50.0, "Pf-α={}", lines[0].1);
    }

    #[test]
    fn series_wavelengths_decrease_with_n() {
        // Within any series, wavelengths should decrease with increasing n
        let lines = balmer_series(1, 10).unwrap();
        for i in 1..lines.len() {
            assert!(
                lines[i].1 < lines[i - 1].1,
                "Wavelength should decrease: {} vs {}",
                lines[i].1,
                lines[i - 1].1
            );
        }
    }

    #[test]
    fn series_invalid_n_max() {
        assert!(lyman_series(1, 1).is_err());
        assert!(balmer_series(1, 2).is_err());
        assert!(paschen_series(1, 3).is_err());
        assert!(brackett_series(1, 4).is_err());
        assert!(pfund_series(1, 5).is_err());
    }

    // --- Dirac energy tests ---

    #[test]
    fn dirac_energy_hydrogen_1s() {
        // H 1s1/2 (n=1, j=1/2 -> two_j=1)
        // Should be very close to electron rest mass (slightly less due to binding)
        let e = dirac_energy_mev(1, 1, 1).unwrap();
        let me = ELECTRON_MASS_MEV;
        assert!(e < me, "Dirac energy should be less than rest mass");
        assert!(e > 0.510_9, "Dirac energy should be close to rest mass");
    }

    #[test]
    fn dirac_binding_hydrogen_1s() {
        // H 1s1/2 binding energy should be ~13.6 eV
        let binding = dirac_binding_energy_ev(1, 1, 1).unwrap();
        assert!(
            (binding - 13.6).abs() < 0.1,
            "H 1s binding={binding} eV, expected ~13.6"
        );
    }

    #[test]
    fn dirac_binding_hydrogen_2s_2p() {
        // 2s1/2 and 2p1/2 should have same Dirac energy (same j=1/2)
        let e_2s = dirac_energy_mev(1, 2, 1).unwrap();
        let e_2p = dirac_energy_mev(1, 2, 1).unwrap(); // same n,j
        assert!(
            (e_2s - e_2p).abs() < 1e-15,
            "2s1/2 and 2p1/2 have same Dirac energy"
        );

        // 2p3/2 should have slightly different energy
        let e_2p32 = dirac_energy_mev(1, 2, 3).unwrap();
        assert!(
            (e_2s - e_2p32).abs() > 1e-12,
            "2s1/2 and 2p3/2 should differ"
        );
    }

    #[test]
    fn dirac_binding_scales_with_z() {
        // Binding energy should scale roughly as Z² for low Z
        let h = dirac_binding_energy_ev(1, 1, 1).unwrap();
        let he = dirac_binding_energy_ev(2, 1, 1).unwrap();
        let ratio = he / h;
        assert!(
            (ratio - 4.0).abs() < 0.01,
            "He/H binding ratio={ratio}, expected ~4"
        );
    }

    #[test]
    fn dirac_energy_invalid_n() {
        assert!(dirac_energy_mev(1, 0, 1).is_err());
    }

    #[test]
    fn dirac_energy_invalid_two_j_even() {
        assert!(dirac_energy_mev(1, 1, 2).is_err());
    }

    #[test]
    fn dirac_energy_invalid_two_j_zero() {
        assert!(dirac_energy_mev(1, 1, 0).is_err());
    }

    #[test]
    fn dirac_energy_high_z() {
        // Uranium Z=92: should still be valid for 1s1/2
        let e = dirac_energy_mev(92, 1, 1).unwrap();
        let me = ELECTRON_MASS_MEV;
        // For high Z, binding is significant
        assert!(e < me, "Bound state energy < rest mass");
        assert!(e > 0.0, "Energy should be positive");
    }

    // --- Relativistic correction tests ---

    #[test]
    fn relativistic_correction_hydrogen_positive() {
        // For hydrogen ground state, relativistic correction should be positive
        // (Dirac binding > non-relativistic binding)
        let corr = relativistic_correction_ev(1, 1, 1).unwrap();
        assert!(
            corr.abs() < 0.01,
            "H 1s correction={corr} eV, should be small"
        );
    }

    #[test]
    fn relativistic_correction_increases_with_z() {
        let corr_h = relativistic_correction_ev(1, 1, 1).unwrap().abs();
        let corr_fe = relativistic_correction_ev(26, 1, 1).unwrap().abs();
        assert!(
            corr_fe > corr_h,
            "Relativistic correction should increase with Z"
        );
    }

    #[test]
    fn relativistic_correction_small_for_hydrogen() {
        // For hydrogen, the correction should be much smaller than 13.6 eV
        let corr = relativistic_correction_ev(1, 1, 1).unwrap();
        assert!(
            corr.abs() < 1.0,
            "H correction={corr} eV, should be << 13.6"
        );
    }

    #[test]
    fn relativistic_correction_invalid() {
        assert!(relativistic_correction_ev(1, 0, 1).is_err());
    }

    // --- Hyperfine splitting tests ---

    #[test]
    fn hyperfine_hydrogen_1s_21cm() {
        // H 1s hyperfine: ~5.88e-6 eV -> 1420 MHz (21 cm line)
        let g_p = crate::constants::PROTON_G_FACTOR;
        let hfs = hyperfine_splitting_ev(1, 1, g_p);
        // The formula gives an approximate value; check order of magnitude
        assert!(
            hfs > 1e-6 && hfs < 1e-5,
            "H 1s HFS={hfs} eV, expected ~5.88e-6"
        );
    }

    #[test]
    fn hyperfine_scales_with_z3() {
        let g_p = crate::constants::PROTON_G_FACTOR;
        let h = hyperfine_splitting_ev(1, 1, g_p);
        let he = hyperfine_splitting_ev(2, 1, g_p);
        let ratio = he / h;
        assert!(
            (ratio - 8.0).abs() < 0.1,
            "HFS Z³ scaling: ratio={ratio}, expected 8"
        );
    }

    #[test]
    fn hyperfine_decreases_with_n() {
        let g_p = crate::constants::PROTON_G_FACTOR;
        let n1 = hyperfine_splitting_ev(1, 1, g_p);
        let n2 = hyperfine_splitting_ev(1, 2, g_p);
        assert!(n1 > n2, "HFS should decrease with n");
    }

    #[test]
    fn hyperfine_zero_n() {
        let hfs = hyperfine_splitting_ev(1, 0, 5.0);
        assert!((hfs).abs() < 1e-30, "n=0 should return 0");
    }

    // --- Electron g-factor tests ---

    #[test]
    fn electron_g_factor_value() {
        let g = electron_g_factor();
        assert!(
            (g - 2.002_319_304).abs() < 1e-6,
            "g_e={g}, expected ~2.00232"
        );
    }

    #[test]
    fn electron_g_factor_greater_than_2() {
        let g = electron_g_factor();
        assert!(g > 2.0, "g_e should be > 2 due to QED corrections");
    }

    #[test]
    fn bound_g_factor_hydrogen() {
        // For Z=1, bound g-factor should be very close to free g-factor
        let g_bound = bound_electron_g_factor(1);
        let g_free = electron_g_factor();
        assert!(
            (g_bound - g_free).abs() < 0.001,
            "H bound g={g_bound}, free g={g_free}"
        );
    }

    #[test]
    fn bound_g_factor_decreases_with_z() {
        let g1 = bound_electron_g_factor(1);
        let g50 = bound_electron_g_factor(50);
        assert!(
            g1 > g50,
            "Bound g-factor should decrease with Z: g(1)={g1}, g(50)={g50}"
        );
    }

    #[test]
    fn bound_g_factor_supercritical() {
        // Z > 137 is supercritical, should return 0
        let g = bound_electron_g_factor(200);
        assert!((g).abs() < 1e-30, "Supercritical should return 0");
    }

    // --- Anomalous Zeeman tests ---

    #[test]
    fn anomalous_zeeman_zero_field() {
        let de = anomalous_zeeman_splitting_ev(0.0, 0, 1);
        assert!((de).abs() < 1e-30, "No splitting in zero field");
    }

    #[test]
    fn anomalous_zeeman_proportional_to_field() {
        let de1 = anomalous_zeeman_splitting_ev(1.0, 0, 1);
        let de2 = anomalous_zeeman_splitting_ev(2.0, 0, 1);
        assert!(
            (de2 - 2.0 * de1).abs() < 1e-15,
            "Anomalous Zeeman should be linear in B"
        );
    }

    #[test]
    fn anomalous_zeeman_differs_from_normal() {
        // The anomalous Zeeman uses g_e != 2, so it should differ
        // from the simple formula with g=2
        let anom = anomalous_zeeman_splitting_ev(1.0, 0, 1);
        let mu_b = crate::constants::BOHR_MAGNETON_EV_T;
        let normal = (0.0 + 2.0 * 0.5) * mu_b * 1.0; // ml=0, ms=+1/2, g=2
        assert!(
            (anom - normal).abs() > 1e-8,
            "Anomalous should differ from normal Zeeman"
        );
    }

    #[test]
    fn anomalous_zeeman_spin_up_vs_down() {
        let up = anomalous_zeeman_splitting_ev(1.0, 0, 1);
        let down = anomalous_zeeman_splitting_ev(1.0, 0, -1);
        assert!(
            (up + down).abs() < 1e-15,
            "Spin-up and spin-down should be symmetric for ml=0"
        );
    }

    // --- Breit interaction tests ---

    #[test]
    fn breit_leading_order_value() {
        // α² Z³ / 4 Hartree for Z = 2: 2 α² E_h.
        let de = breit_interaction_ev(2);
        let expected = 2.0 * FINE_STRUCTURE * FINE_STRUCTURE * 2.0 * RYDBERG_EV;
        assert!((de - expected).abs() < 1e-15, "He Breit={de}");
        assert!(de > 0.0);
    }

    #[test]
    fn breit_scales_with_z3() {
        let ratio = breit_interaction_ev(20) / breit_interaction_ev(10);
        assert!((ratio - 8.0).abs() < 1e-12, "ratio={ratio}");
    }

    #[test]
    fn breit_below_total_two_electron_energy_uranium() {
        // Measured total two-electron contribution for He-like U: 2248 ± 9 eV
        // (Gumberidze et al., PRL 92, 203004 (2004)).
        assert!(breit_interaction_ev(92) < 2248.0);
    }

    // --- Reference-value tests ---

    fn mhz(ev: f64) -> f64 {
        ev / H_EV_S / 1e6
    }

    #[test]
    fn radial_wavefunctions_normalized() {
        for n in 1..=10_u32 {
            for l in 0..n {
                let steps = 40_000;
                let r_max = 4.0 * f64::from(n * n) + 40.0;
                let h = r_max / f64::from(steps);
                let f = |r: f64| radial_probability_density(1, n, l, r).unwrap();
                let mut acc = f(0.0) + f(r_max);
                for k in 1..steps {
                    acc += if k % 2 == 1 { 4.0 } else { 2.0 } * f(f64::from(k) * h);
                }
                let norm = acc * h / 3.0;
                assert!((norm - 1.0).abs() < 1e-9, "n={n} l={l}: norm={norm}");
            }
        }
    }

    #[test]
    fn radial_r31_matches_textbook() {
        // R_31 = (8/(27√6)) (1 − r/6) r e^{−r/3} for Z = 1 (Griffiths eq. 4.89).
        for r in [0.5, 2.0, 6.0, 9.0] {
            let expected =
                8.0 / (27.0 * libm::sqrt(6.0)) * (1.0 - r / 6.0) * r * libm::exp(-r / 3.0);
            let got = radial_wavefunction(1, 3, 1, r).unwrap();
            assert!((got - expected).abs() < 1e-15, "r={r}: {got} vs {expected}");
        }
    }

    #[test]
    fn einstein_a_matches_nist_with_reduced_mass() {
        // NIST ASD hydrogen A-values (s⁻¹), which include the reduced mass.
        let mu = reduced_mass_factor(1, 1).unwrap();
        for (nu, lu, nl, ll, nist) in [
            (2, 1, 1, 0, 6.2649e8),
            (3, 1, 1, 0, 1.6725e8),
            (3, 1, 2, 0, 2.2448e7),
            (3, 0, 2, 1, 6.3143e6),
            (3, 2, 2, 1, 6.4651e7),
        ] {
            let a = einstein_a_coefficient(1, nu, lu, nl, ll).unwrap() * mu;
            assert!(
                (a - nist).abs() / nist < 2e-4,
                "{nu}{lu}->{nl}{ll}: {a} vs {nist}"
            );
        }
    }

    #[test]
    fn einstein_a_scales_as_z4() {
        let h = einstein_a_coefficient(1, 2, 1, 1, 0).unwrap();
        let he = einstein_a_coefficient(2, 2, 1, 1, 0).unwrap();
        assert!((he / h - 16.0).abs() < 1e-9, "ratio={}", he / h);
    }

    #[test]
    fn vacuum_polarization_hydrogen_2s() {
        // H 2S: leading Uehling term −27.13 MHz plus α(Zα)⁵ term +0.24 MHz.
        let vp = mhz(vacuum_polarization_ev(1, 2, 0));
        assert!((vp + 26.886).abs() < 0.01, "VP={vp} MHz");
    }

    #[test]
    fn lamb_shift_hydrogen_classic() {
        // 2S½ − 2P½ measured 1057.845 MHz; 1S Lamb shift 8172.9 MHz.
        let split =
            mhz(lamb_shift_nlj_ev(1, 2, 0, 1).unwrap() - lamb_shift_nlj_ev(1, 2, 1, 1).unwrap());
        assert!(
            (split - 1057.845).abs() / 1057.845 < 3e-3,
            "2S-2P={split} MHz"
        );
        let one_s = mhz(lamb_shift_nlj_ev(1, 1, 0, 1).unwrap());
        assert!((one_s - 8172.9).abs() / 8172.9 < 3e-3, "1S={one_s} MHz");
        let p32 = mhz(lamb_shift_nlj_ev(1, 2, 1, 3).unwrap());
        assert!((p32 - 12.8).abs() < 0.5, "2P3/2={p32} MHz");
    }

    #[test]
    fn lamb_shift_helium_ion() {
        // He⁺ 2S½ − 2P½ measured 14041.13 MHz.
        let split =
            mhz(lamb_shift_nlj_ev(2, 2, 0, 1).unwrap() - lamb_shift_nlj_ev(2, 2, 1, 1).unwrap());
        assert!(
            (split - 14_041.13).abs() / 14_041.13 < 0.04,
            "He+ 2S-2P={split} MHz"
        );
    }

    #[test]
    fn lamb_shift_rejects_untabulated() {
        assert!(lamb_shift_nlj_ev(1, 5, 0, 1).is_err());
        assert!(lamb_shift_nlj_ev(1, 2, 1, 5).is_err());
    }

    #[test]
    fn hyperfine_deuterium_includes_spin_factor() {
        let hfs = mhz(hyperfine_splitting_spin_ev(
            1,
            1,
            crate::constants::DEUTERON_G_FACTOR,
            1.0,
        ));
        assert!((hfs - 327.384).abs() / 327.384 < 1e-3, "D HFS={hfs} MHz");
        let h = mhz(hyperfine_splitting_ev(
            1,
            1,
            crate::constants::PROTON_G_FACTOR,
        ));
        assert!((h - 1420.406).abs() / 1420.406 < 1e-3, "H HFS={h} MHz");
    }

    #[test]
    fn h_alpha_vacuum_with_reduced_mass() {
        let lambda = spectral_line_vacuum_nm(1, 1, 2, 3).unwrap();
        assert!((lambda - 656.4696).abs() < 1e-3, "Hα={lambda} nm");
    }

    #[test]
    fn fine_structure_uses_codata_rydberg() {
        let e = hydrogen_level_energy_ev(1, 1, 1).unwrap();
        let dirac = -dirac_binding_energy_ev(1, 1, 1).unwrap();
        assert!((e - dirac).abs() < 1e-8, "first-order {e} vs Dirac {dirac}");
        assert!(hydrogen_level_energy_ev(1, 2, 2).is_err());
    }

    #[test]
    fn dirac_rejects_j_above_n() {
        assert!(dirac_energy_mev(1, 1, 3).is_err());
    }

    #[test]
    fn ionization_energy_reference_values() {
        // NIST ASD.
        assert!((ionization_energy_ev(43).unwrap() - 7.119_38).abs() < 1e-9);
        assert!((ionization_energy_ev(96).unwrap() - 5.992_241).abs() < 1e-9);
        // Smits et al. 2023.
        assert!((ionization_energy_ev(112).unwrap() - 12.02).abs() < 1e-9);
    }

    #[test]
    fn electron_affinity_reference_values() {
        assert_eq!(
            electron_affinity(90).unwrap(),
            ElectronAffinity::Bound(0.607_69)
        );
        assert_eq!(
            electron_affinity(92).unwrap(),
            ElectronAffinity::Bound(0.314_97)
        );
        assert_eq!(electron_affinity(2).unwrap(), ElectronAffinity::Unbound);
        assert_eq!(electron_affinity(104).unwrap(), ElectronAffinity::Unknown);
        assert!((electron_affinity_ev(59).unwrap() - 0.109_23).abs() < 1e-12);
    }

    #[test]
    fn superheavy_configurations_follow_madelung() {
        let ds = format_configuration_short(&electron_configuration(110).unwrap(), 110);
        assert_eq!(ds, "[Rn] 7s2 5f14 6d8");
        let rg = format_configuration_short(&electron_configuration(111).unwrap(), 111);
        assert_eq!(rg, "[Rn] 7s2 5f14 6d9");
    }

    #[test]
    fn stark_and_lande_reject_impossible_states() {
        assert!(stark_shift_hydrogen_ev(2, 5, 1e6).abs() < 1e-30);
        assert!(lande_g_factor(0, 3).abs() < 1e-30);
    }
}
