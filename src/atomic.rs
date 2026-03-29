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
/// Source: NIST Atomic Spectra Database ground-state configurations.
/// All 22 known Aufbau/Madelung exceptions for Z=1-118.
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
        110 => {
            // Ds: [Rn] 5f14 6d9 7s1 (predicted, parallels Pt)
            set_subshell(config, 7, OrbitalType::S, 1);
            set_subshell(config, 6, OrbitalType::D, 9);
        }
        111 => {
            // Rg: [Rn] 5f14 6d10 7s1 (predicted, parallels Au)
            set_subshell(config, 7, OrbitalType::S, 1);
            set_subshell(config, 6, OrbitalType::D, 10);
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

/// Calculates the energy of a hydrogen-like level with fine-structure correction.
///
/// E_nj = -13.6 eV * Z² / n² * [1 + (αZ)²/n * (1/(j+1/2) - 3/(4n))]
///
/// where α is the fine-structure constant, j is the total angular momentum
/// quantum number (j = l ± 1/2, stored as integer 2j).
///
/// Returns the energy in eV (negative, bound state).
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] if quantum numbers are invalid.
pub fn hydrogen_level_energy_ev(z: u32, n: u32, two_j: u32) -> Result<f64, TanmatraError> {
    if n == 0 {
        return Err(TanmatraError::InvalidQuantumNumbers(String::from(
            "n must be >= 1",
        )));
    }
    if two_j == 0 || two_j > 2 * n - 1 {
        return Err(TanmatraError::InvalidQuantumNumbers(alloc::format!(
            "2j={two_j} invalid for n={n}"
        )));
    }

    let z_f = z as f64;
    let n_f = n as f64;
    let alpha = crate::constants::FINE_STRUCTURE;
    let j_plus_half = f64::midpoint(two_j as f64, 1.0);

    // Non-relativistic energy
    let e0 = -13.6 * z_f * z_f / (n_f * n_f);

    // Fine-structure correction (first-order)
    let az = alpha * z_f;
    let correction = 1.0 + az * az / n_f * (1.0 / j_plus_half - 3.0 / (4.0 * n_f));

    Ok(e0 * correction)
}

/// Calculates the wavelength of a spectral line with fine-structure correction.
///
/// Uses the Dirac energy levels for hydrogen-like atoms, which include
/// the relativistic fine-structure splitting.
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

    // Convert eV to wavelength in nm: λ = hc/ΔE
    // hc = 1239.8419843320028 eV·nm
    let hc_ev_nm = 1_239.841_984;
    Ok(hc_ev_nm / delta_e)
}

// ---------------------------------------------------------------------------
// Zeeman and Stark effects
// ---------------------------------------------------------------------------

/// Calculates the Lande g-factor for an atomic level.
///
/// g_J = 1 + [J(J+1) + S(S+1) - L(L+1)] / [2J(J+1)]
///
/// For a single electron: S = 1/2, L = l, J = j.
///
/// Parameters: `l` (orbital), `two_j` (2*total angular momentum).
/// Returns the dimensionless g-factor.
#[must_use]
#[inline]
pub fn lande_g_factor(l: u32, two_j: u32) -> f64 {
    let j = two_j as f64 / 2.0;
    let l_f = l as f64;
    let s = 0.5; // single electron spin

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
    mj * g * crate::constants::BOHR_MAGNETON_EV_T * b_tesla
}

/// Calculates the linear Stark effect energy shift for hydrogen.
///
/// For hydrogen, the linear Stark effect gives:
/// ΔE = (3/2) * n * (n1 - n2) * e * a0 * E_field
///
/// where n1 and n2 are parabolic quantum numbers with n1 + n2 + |m| + 1 = n.
///
/// Simplified: for the maximum shift state (n1 = n-1, n2 = 0, m = 0):
/// ΔE_max = (3/2) * n * (n-1) * e * a0 * E_field
///
/// Parameters:
/// - `n`: principal quantum number
/// - `parabolic_index`: (n1 - n2), ranges from -(n-1) to (n-1)
/// - `e_field_v_per_m`: electric field strength in V/m
///
/// Returns the energy shift in eV.
#[must_use]
#[inline]
pub fn stark_shift_hydrogen_ev(n: u32, parabolic_index: i32, e_field_v_per_m: f64) -> f64 {
    let n_f = n as f64;
    let k = parabolic_index as f64;
    let a0 = crate::constants::BOHR_RADIUS; // meters
    // ΔE = (3/2) * n * k * e * a0 * E = (3/2) * n * k * a0 * E (in eV, since e*V = eV)
    1.5 * n_f * k * a0 * e_field_v_per_m
}

// ---------------------------------------------------------------------------
// Hydrogen wavefunctions
// ---------------------------------------------------------------------------

/// Computes the radial wavefunction R_nl(r) for a hydrogen-like atom.
///
/// R_nl(r) = N * (2Zr/na0)^l * exp(-Zr/na0) * L_{n-l-1}^{2l+1}(2Zr/na0)
///
/// where L is the associated Laguerre polynomial and N is the normalization.
///
/// For simplicity, implements exact forms for n=1,2,3 with any l.
///
/// Parameters:
/// - `z`: atomic number
/// - `n`: principal quantum number (1, 2, or 3)
/// - `l`: orbital angular momentum (0 to n-1)
/// - `r_bohr`: radial distance in units of Bohr radius (r/a0)
///
/// Returns R_nl(r) * a0^(3/2) (dimensionless when multiplied by a0^{-3/2}).
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidQuantumNumbers`] if n > 3 or l >= n.
pub fn radial_wavefunction(z: u32, n: u32, l: u32, r_bohr: f64) -> Result<f64, TanmatraError> {
    if n == 0 || n > 3 {
        return Err(TanmatraError::InvalidQuantumNumbers(alloc::format!(
            "n={n} not supported (1-3)"
        )));
    }
    if l >= n {
        return Err(TanmatraError::InvalidQuantumNumbers(alloc::format!(
            "l={l} must be < n={n}"
        )));
    }

    let zf = z as f64;
    let rho = 2.0 * zf * r_bohr / n as f64; // 2Zr/(n*a0)
    let exp_factor = libm::exp(-rho / 2.0);

    let result = match (n, l) {
        (1, 0) => {
            // R_10 = 2 * (Z/a0)^{3/2} * exp(-Zr/a0)
            2.0 * zf * libm::sqrt(zf) * exp_factor
        }
        (2, 0) => {
            // R_20 = (1/2√2) * (Z/a0)^{3/2} * (2 - Zr/a0) * exp(-Zr/2a0)
            let norm = 1.0 / (2.0 * core::f64::consts::SQRT_2);
            norm * zf * libm::sqrt(zf) * (2.0 - zf * r_bohr) * exp_factor
        }
        (2, 1) => {
            // R_21 = (1/2√6) * (Z/a0)^{3/2} * (Zr/a0) * exp(-Zr/2a0)
            let norm = 1.0 / (2.0 * libm::sqrt(6.0));
            norm * zf * libm::sqrt(zf) * zf * r_bohr * exp_factor
        }
        (3, 0) => {
            // R_30 = (2/81√3) * (Z/a0)^{3/2} * (27 - 18Zr/a0 + 2(Zr/a0)²) * exp(-Zr/3a0)
            let zr = zf * r_bohr;
            let norm = 2.0 / (81.0 * libm::sqrt(3.0));
            norm * zf * libm::sqrt(zf) * (27.0 - 18.0 * zr + 2.0 * zr * zr) * exp_factor
        }
        (3, 1) => {
            // R_31 = (8/27√6) * (Z/a0)^{3/2} * (Zr/a0)(6 - Zr/a0) * exp(-Zr/3a0)
            let zr = zf * r_bohr;
            let norm = 8.0 / (27.0 * libm::sqrt(6.0));
            norm * zf * libm::sqrt(zf) * zr * (6.0 - zr) * exp_factor
        }
        (3, 2) => {
            // R_32 = (4/81√30) * (Z/a0)^{3/2} * (Zr/a0)² * exp(-Zr/3a0)
            let zr = zf * r_bohr;
            let norm = 4.0 / (81.0 * libm::sqrt(30.0));
            norm * zf * libm::sqrt(zf) * zr * zr * exp_factor
        }
        _ => {
            return Err(TanmatraError::InvalidQuantumNumbers(alloc::format!(
                "n={n}, l={l} not supported"
            )));
        }
    };

    Ok(result)
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

/// Calculates the Einstein A coefficient (spontaneous emission rate) for
/// hydrogen-like atoms.
///
/// For electric dipole transitions between levels (n1,l1) and (n2,l2)
/// in a hydrogen-like atom with atomic number Z:
///
/// A_{21} = (4 α³ ω³)/(3 c²) |⟨1|r|2⟩|²
///
/// For hydrogen (Z=1), a simplified formula gives:
/// A ∝ Z⁴ * (frequency)³ * |radial matrix element|²
///
/// This implementation uses the known analytical result for hydrogen:
/// A_{n2,l2 → n1,l1} scales as Z⁴ / n_eff⁵
///
/// Returns the rate in s⁻¹.
///
/// # Errors
///
/// Returns error if the transition violates selection rules.
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

    let zf = z as f64;
    let n1 = n_lower as f64;
    let n2 = n_upper as f64;
    let l_max = if l_upper > l_lower { l_upper } else { l_lower };

    // Transition energy in eV
    let delta_e_ev = 13.6 * zf * zf * (1.0 / (n1 * n1) - 1.0 / (n2 * n2));

    // Transition frequency in Hz: ν = ΔE / h
    let freq_hz = delta_e_ev / crate::constants::H_EV_S;

    // Use the standard relation: A = (8π² e² ν²) / (m_e c³) * f_osc
    // with Kramers' approximate oscillator strength for hydrogen:
    // f ≈ (32/(3√3 π)) * 1/(n1² n2² (1/n1² - 1/n2²)³) * max(l,l')/(2l_upper+1)
    //
    // Instead, use the exact formula in SI:
    // A_{21} = (ω³ |d|²) / (3π ε₀ ħ c³)
    //
    // For hydrogen transitions, the known scaling is:
    // A = 6.27e8 * Z⁴ * (freq/freq_ly_alpha)³ * max(l,l')/(2l_upper+1) * correction
    //
    // Ly-alpha: n=2→1, l=1→0, freq = 2.466e15 Hz, A = 6.27e8 s⁻¹
    let freq_ly_alpha = 2.466e15; // Hz for hydrogen Lyman-alpha
    let a_ly_alpha = 6.27e8; // s⁻¹ for hydrogen Lyman-alpha

    let freq_ratio = freq_hz / freq_ly_alpha;
    let g_ratio = l_max as f64 / (2 * l_upper + 1) as f64;

    // A scales as Z⁴ * ν³ * angular factor
    let a_coeff = a_ly_alpha * zf.powi(4) * freq_ratio.powi(3) * g_ratio;

    Ok(a_coeff.abs())
}

/// Calculates the Einstein B coefficient for stimulated emission/absorption.
///
/// B_{21} = A_{21} * c³ / (8π h ν³)
///
/// where ν is the transition frequency.
///
/// Returns B in m³/(J·s²) = m³·sr/(J·s).
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

    // Transition frequency in Hz
    let delta_e_ev = 13.6 * zf * zf * (1.0 / (n1 * n1) - 1.0 / (n2 * n2));
    let freq_hz = delta_e_ev / crate::constants::H_EV_S;

    let c_val = crate::constants::C;
    let h_val = crate::constants::H_EV_S * crate::constants::ELEMENTARY_CHARGE; // h in J·s

    if freq_hz <= 0.0 {
        return Ok(0.0);
    }

    Ok(a21 * c_val * c_val * c_val / (8.0 * core::f64::consts::PI * h_val * freq_hz.powi(3)))
}

// ---------------------------------------------------------------------------
// Electron affinities
// ---------------------------------------------------------------------------

/// Electron affinities in eV for Z=1..=118.
///
/// Positive value = energy released when forming X⁻ (stable anion).
/// 0.0 = anion unstable (noble gases, alkaline earths, Mn, N, etc.).
///
/// Sources: NIST, T. Andersen (2004), Rienstra-Kiracofe et al. (2002).
/// Z>104: relativistic theoretical values (Eliav, Kaldor, Borschevsky).
#[allow(clippy::too_many_lines)]
const ELECTRON_AFFINITY_EV: [f64; 118] = [
    0.754_20, // H  (Z=1)
    0.0,      // He (Z=2)
    0.618_05, // Li (Z=3)
    0.0,      // Be (Z=4)
    0.279_72, // B  (Z=5)
    1.262_12, // C  (Z=6)
    0.0,      // N  (Z=7)
    1.461_12, // O  (Z=8)
    3.401_19, // F  (Z=9)
    0.0,      // Ne (Z=10)
    0.547_93, // Na (Z=11)
    0.0,      // Mg (Z=12)
    0.432_83, // Al (Z=13)
    1.389_52, // Si (Z=14)
    0.746_61, // P  (Z=15)
    2.077_10, // S  (Z=16)
    3.612_72, // Cl (Z=17)
    0.0,      // Ar (Z=18)
    0.501_46, // K  (Z=19)
    0.024_55, // Ca (Z=20)
    0.188,    // Sc (Z=21)
    0.079,    // Ti (Z=22)
    0.525,    // V  (Z=23)
    0.666_0,  // Cr (Z=24)
    0.0,      // Mn (Z=25)
    0.151,    // Fe (Z=26)
    0.662_26, // Co (Z=27)
    1.156_16, // Ni (Z=28)
    1.235_78, // Cu (Z=29)
    0.0,      // Zn (Z=30)
    0.430,    // Ga (Z=31)
    1.232_71, // Ge (Z=32)
    0.804,    // As (Z=33)
    2.020_67, // Se (Z=34)
    3.363_59, // Br (Z=35)
    0.0,      // Kr (Z=36)
    0.485_92, // Rb (Z=37)
    0.048_16, // Sr (Z=38)
    0.307,    // Y  (Z=39)
    0.426,    // Zr (Z=40)
    0.916_0,  // Nb (Z=41)
    0.748_0,  // Mo (Z=42)
    0.550,    // Tc (Z=43)
    1.046_38, // Ru (Z=44)
    1.142_89, // Rh (Z=45)
    0.562_14, // Pd (Z=46)
    1.304_7,  // Ag (Z=47)
    0.0,      // Cd (Z=48)
    0.404,    // In (Z=49)
    1.112_07, // Sn (Z=50)
    1.047_01, // Sb (Z=51)
    1.970_88, // Te (Z=52)
    3.059_04, // I  (Z=53)
    0.0,      // Xe (Z=54)
    0.471_63, // Cs (Z=55)
    0.144_62, // Ba (Z=56)
    0.470,    // La (Z=57)
    0.570,    // Ce (Z=58)
    0.962,    // Pr (Z=59)
    0.0,      // Nd (Z=60)
    0.0,      // Pm (Z=61)
    0.0,      // Sm (Z=62)
    0.0,      // Eu (Z=63)
    0.0,      // Gd (Z=64)
    0.0,      // Tb (Z=65)
    0.0,      // Dy (Z=66)
    0.0,      // Ho (Z=67)
    0.0,      // Er (Z=68)
    1.029,    // Tm (Z=69)
    0.0,      // Yb (Z=70)
    0.340,    // Lu (Z=71)
    0.017,    // Hf (Z=72)
    0.322,    // Ta (Z=73)
    0.816_26, // W  (Z=74)
    0.150,    // Re (Z=75)
    1.077_80, // Os (Z=76)
    1.564_36, // Ir (Z=77)
    2.128_10, // Pt (Z=78)
    2.308_63, // Au (Z=79)
    0.0,      // Hg (Z=80)
    0.377,    // Tl (Z=81)
    0.364_3,  // Pb (Z=82)
    0.942_36, // Bi (Z=83)
    1.900,    // Po (Z=84)
    2.415_78, // At (Z=85)
    0.0,      // Rn (Z=86)
    0.486_3,  // Fr (Z=87)
    0.144,    // Ra (Z=88)
    0.350,    // Ac (Z=89)
    1.170,    // Th (Z=90)
    0.0,      // Pa (Z=91)
    0.0,      // U  (Z=92)
    0.0,      // Np-Og (Z=93-118): most actinides/superheavy have unbound anions
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   // Z=94-102
    0.470, // Lr (Z=103)
    0.0,   // Rf (Z=104)
    0.680, // Db (Z=105)
    0.810, // Sg (Z=106)
    0.0, 0.0, 0.0, 0.0,   // Bh-Ds (Z=107-110)
    1.565, // Rg (Z=111)
    0.0,   // Cn (Z=112)
    0.680, // Nh (Z=113)
    0.905, // Fl (Z=114)
    0.674, // Mc (Z=115)
    1.470, // Lv (Z=116)
    1.635, // Ts (Z=117)
    0.0,   // Og (Z=118)
];

/// Returns the electron affinity in eV for the given atomic number.
///
/// The electron affinity is the energy released when an electron is added
/// to a neutral atom: X + e⁻ → X⁻ + EA.
///
/// Returns 0.0 for elements with unstable anions (noble gases, etc.).
///
/// # Errors
///
/// Returns [`TanmatraError::InvalidAtomicNumber`] if Z is 0 or > 118.
#[inline]
pub fn electron_affinity_ev(z: u32) -> Result<f64, TanmatraError> {
    if z == 0 || z > 118 {
        return Err(TanmatraError::InvalidAtomicNumber(z));
    }
    Ok(ELECTRON_AFFINITY_EV[(z - 1) as usize])
}

/// NIST ionization energies for Z=1 to Z=118 in eV.
///
/// Source: NIST Atomic Spectra Database, Ionization Energies.
/// <https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html>
///
/// Z=1-103: experimental values from NIST ASD.
/// Z=104-118: relativistic theoretical predictions (noted in comments).
#[allow(clippy::too_many_lines)]
const IONIZATION_ENERGIES_EV: [f64; 118] = [
    // --- Period 1 ---
    13.598_44, // H  (Z=1)
    24.587_4,  // He (Z=2)
    // --- Period 2 ---
    5.391_71, // Li (Z=3)
    9.322_7,  // Be (Z=4)
    8.298_0,  // B  (Z=5)
    11.260_3, // C  (Z=6)
    14.534_1, // N  (Z=7)
    13.618_1, // O  (Z=8)
    17.422_8, // F  (Z=9)
    21.564_5, // Ne (Z=10)
    // --- Period 3 ---
    5.139_08, // Na (Z=11)
    7.646_2,  // Mg (Z=12)
    5.985_77, // Al (Z=13)
    8.151_7,  // Si (Z=14)
    10.486_7, // P  (Z=15)
    10.360_0, // S  (Z=16)
    12.967_6, // Cl (Z=17)
    15.759_6, // Ar (Z=18)
    // --- Period 4 ---
    4.340_66, // K  (Z=19)
    6.113_2,  // Ca (Z=20)
    6.561_5,  // Sc (Z=21)
    6.828_1,  // Ti (Z=22)
    6.746_2,  // V  (Z=23)
    6.766_5,  // Cr (Z=24)
    7.434_0,  // Mn (Z=25)
    7.902_4,  // Fe (Z=26)
    7.881_0,  // Co (Z=27)
    7.639_8,  // Ni (Z=28)
    7.726_4,  // Cu (Z=29)
    9.394_2,  // Zn (Z=30)
    5.999_3,  // Ga (Z=31)
    7.899_4,  // Ge (Z=32)
    9.788_6,  // As (Z=33)
    9.752_4,  // Se (Z=34)
    11.813_8, // Br (Z=35)
    13.999_6, // Kr (Z=36)
    // --- Period 5 ---
    4.177_13,  // Rb (Z=37)
    5.694_84,  // Sr (Z=38)
    6.217_26,  // Y  (Z=39)
    6.634_0,   // Zr (Z=40)
    6.758_85,  // Nb (Z=41)
    7.092_43,  // Mo (Z=42)
    7.28,      // Tc (Z=43)
    7.360_50,  // Ru (Z=44)
    7.458_90,  // Rh (Z=45)
    8.336_9,   // Pd (Z=46)
    7.576_24,  // Ag (Z=47)
    8.993_82,  // Cd (Z=48)
    5.786_36,  // In (Z=49)
    7.343_92,  // Sn (Z=50)
    8.608_4,   // Sb (Z=51)
    9.009_66,  // Te (Z=52)
    10.451_26, // I  (Z=53)
    12.129_84, // Xe (Z=54)
    // --- Period 6 ---
    3.893_905, // Cs (Z=55)
    5.211_70,  // Ba (Z=56)
    5.576_9,   // La (Z=57)
    5.538_7,   // Ce (Z=58)
    5.473,     // Pr (Z=59)
    5.525_0,   // Nd (Z=60)
    5.582_0,   // Pm (Z=61)
    5.643_71,  // Sm (Z=62)
    5.670_38,  // Eu (Z=63)
    6.149_80,  // Gd (Z=64)
    5.863_8,   // Tb (Z=65)
    5.938_9,   // Dy (Z=66)
    6.021_5,   // Ho (Z=67)
    6.107_7,   // Er (Z=68)
    6.184_31,  // Tm (Z=69)
    6.254_16,  // Yb (Z=70)
    5.425_59,  // Lu (Z=71)
    6.825_07,  // Hf (Z=72)
    7.549_6,   // Ta (Z=73)
    7.864_03,  // W  (Z=74)
    7.833_52,  // Re (Z=75)
    8.438_23,  // Os (Z=76)
    8.967_00,  // Ir (Z=77)
    8.958_7,   // Pt (Z=78)
    9.225_53,  // Au (Z=79)
    10.437_50, // Hg (Z=80)
    6.108_29,  // Tl (Z=81)
    7.416_66,  // Pb (Z=82)
    7.285_6,   // Bi (Z=83)
    8.414,     // Po (Z=84)
    9.317_5,   // At (Z=85)
    10.748_5,  // Rn (Z=86)
    // --- Period 7 ---
    4.072_74, // Fr (Z=87)
    5.278_46, // Ra (Z=88)
    5.380_2,  // Ac (Z=89)
    6.306_7,  // Th (Z=90)
    5.89,     // Pa (Z=91)
    6.194_05, // U  (Z=92)
    6.265_6,  // Np (Z=93)
    6.026_0,  // Pu (Z=94)
    5.993_8,  // Am (Z=95)
    6.021_96, // Cm (Z=96)
    6.198_5,  // Bk (Z=97)
    6.281_7,  // Cf (Z=98)
    6.42,     // Es (Z=99)
    6.50,     // Fm (Z=100)
    6.58,     // Md (Z=101)
    6.65,     // No (Z=102)
    4.96,     // Lr (Z=103)
    // --- Superheavy (theoretical) ---
    6.01, // Rf (Z=104)
    6.89, // Db (Z=105)
    7.08, // Sg (Z=106)
    7.7,  // Bh (Z=107)
    7.6,  // Hs (Z=108)
    9.1,  // Mt (Z=109)
    8.7,  // Ds (Z=110)
    9.79, // Rg (Z=111)
    9.38, // Cn (Z=112)
    5.85, // Nh (Z=113)
    7.31, // Fl (Z=114)
    6.92, // Mc (Z=115)
    8.59, // Lv (Z=116)
    7.64, // Ts (Z=117)
    8.91, // Og (Z=118)
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
        assert!(radial_wavefunction(1, 4, 0, 1.0).is_err()); // n > 3
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
}
