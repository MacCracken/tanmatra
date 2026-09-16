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

/// Generates a decay chain from the given nucleus until a nucleus with no
/// known decay is reached or `max_steps` is exceeded.
///
/// Uses the known isotope database to determine decay modes and daughters,
/// following the dominant branch at each step. When a nucleus appears in the
/// database both as an isomer and as a ground state, a freshly produced
/// daughter takes the first listed entry (e.g. Th-234 β⁻ feeds Pa-234m), and
/// after an isomeric transition only the ground state is considered, so the
/// chain continues from — or ends at — the ground state instead of repeating
/// the transition. Entries with no observed decay (infinite half-life) end
/// the chain.
///
/// Returns a list of (Nucleus, DecayMode) pairs for each step.
#[must_use]
pub fn decay_chain(start: &Nucleus, max_steps: usize) -> Vec<(Nucleus, DecayMode)> {
    let mut chain = Vec::new();
    let mut current = *start;
    let mut in_ground_state = false;
    let all = known_isotopes();

    for _ in 0..max_steps {
        let isotope = all
            .iter()
            .find(|iso| iso.nucleus == current && !(in_ground_state && iso.is_isomer));

        let Some(iso) = isotope else { break };
        if !iso.half_life_seconds.is_finite() {
            break;
        }
        let mode = iso.primary_decay;
        let daughter = match mode {
            DecayMode::Alpha => alpha_decay(&current),
            DecayMode::BetaMinus => beta_minus_decay(&current),
            DecayMode::BetaPlus | DecayMode::ElectronCapture => beta_plus_decay(&current),
            DecayMode::IsomericTransition | DecayMode::Gamma => {
                if iso.is_isomer {
                    Ok(current)
                } else {
                    break;
                }
            }
            _ => break, // Fission, proton/neutron emission -- stop chain
        };

        match daughter {
            Ok(d) => {
                chain.push((current, mode));
                in_ground_state = matches!(mode, DecayMode::IsomericTransition | DecayMode::Gamma);
                current = d;
            }
            Err(_) => break,
        }
    }

    chain
}

/// Returns a collection of known radioactive isotopes with real half-lives.
///
/// Half-lives, excitation energies and dominant decay modes are from
/// NUBASE2020 (Kondev et al., Chin. Phys. C 45, 030001 (2021)), with
/// 1 y = 365.2422 d. `primary_decay` is the branch with the largest intensity.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn known_isotopes() -> Vec<Isotope> {
    // NUBASE2020 convention: 1 y = 365.2422 d = 31 556 926 s.
    let seconds_per_year = 365.2422 * 24.0 * 3600.0;
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
        gs(6, 14, "C-14", 5.70e3 * yr, DecayMode::BetaMinus),
        gs(9, 18, "F-18", 109.734 * mn, DecayMode::BetaPlus),
        gs(11, 22, "Na-22", 2.6019 * yr, DecayMode::BetaPlus),
        gs(13, 26, "Al-26", 717.0e3 * yr, DecayMode::BetaPlus),
        gs(14, 32, "Si-32", 157.0 * yr, DecayMode::BetaMinus),
        gs(15, 32, "P-32", 14.269 * dy, DecayMode::BetaMinus),
        gs(16, 35, "S-35", 87.37 * dy, DecayMode::BetaMinus),
        gs(17, 36, "Cl-36", 301.3e3 * yr, DecayMode::BetaMinus),
        gs(18, 39, "Ar-39", 268.0 * yr, DecayMode::BetaMinus),
        gs(19, 40, "K-40", 1.248e9 * yr, DecayMode::BetaMinus),
        gs(20, 41, "Ca-41", 99.4e3 * yr, DecayMode::ElectronCapture),
        gs(24, 51, "Cr-51", 27.7015 * dy, DecayMode::ElectronCapture),
        gs(25, 53, "Mn-53", 3.7e6 * yr, DecayMode::ElectronCapture),
        gs(26, 55, "Fe-55", 2.7562 * yr, DecayMode::ElectronCapture),
        gs(26, 59, "Fe-59", 44.500 * dy, DecayMode::BetaMinus),
        gs(27, 57, "Co-57", 271.811 * dy, DecayMode::ElectronCapture),
        gs(27, 60, "Co-60", 5.2714 * yr, DecayMode::BetaMinus),
        gs(28, 63, "Ni-63", 101.2 * yr, DecayMode::BetaMinus),
        gs(29, 64, "Cu-64", 12.7004 * hr, DecayMode::ElectronCapture), // ε 61.52% (EC 44.0%, β⁺ 17.5%), β⁻ 38.48%
        gs(30, 65, "Zn-65", 243.94 * dy, DecayMode::ElectronCapture),
        gs(31, 67, "Ga-67", 3.2617 * dy, DecayMode::ElectronCapture),
        gs(31, 68, "Ga-68", 67.842 * mn, DecayMode::BetaPlus),
        gs(36, 85, "Kr-85", 10.728 * yr, DecayMode::BetaMinus),
        gs(37, 87, "Rb-87", 49.7e9 * yr, DecayMode::BetaMinus),
        gs(38, 89, "Sr-89", 50.563 * dy, DecayMode::BetaMinus),
        gs(38, 90, "Sr-90", 28.91 * yr, DecayMode::BetaMinus),
        gs(39, 90, "Y-90", 64.05 * hr, DecayMode::BetaMinus),
        gs(40, 95, "Zr-95", 64.032 * dy, DecayMode::BetaMinus),
        gs(41, 95, "Nb-95", 34.991 * dy, DecayMode::BetaMinus),
        gs(42, 99, "Mo-99", 65.932 * hr, DecayMode::BetaMinus),
        gs(43, 99, "Tc-99", 211.1e3 * yr, DecayMode::BetaMinus),
        gs(44, 103, "Ru-103", 39.245 * dy, DecayMode::BetaMinus),
        gs(44, 106, "Ru-106", 371.8 * dy, DecayMode::BetaMinus),
        gs(48, 109, "Cd-109", 461.3 * dy, DecayMode::ElectronCapture),
        gs(49, 111, "In-111", 2.8048 * dy, DecayMode::ElectronCapture),
        gs(50, 113, "Sn-113", 115.08 * dy, DecayMode::ElectronCapture),
        gs(51, 125, "Sb-125", 2.7576 * yr, DecayMode::BetaMinus),
        gs(52, 132, "Te-132", 3.204 * dy, DecayMode::BetaMinus),
        gs(53, 123, "I-123", 13.2232 * hr, DecayMode::ElectronCapture),
        gs(53, 125, "I-125", 59.392 * dy, DecayMode::ElectronCapture),
        gs(53, 129, "I-129", 16.14e6 * yr, DecayMode::BetaMinus),
        gs(53, 131, "I-131", 8.0249 * dy, DecayMode::BetaMinus),
        gs(55, 134, "Cs-134", 2.0650 * yr, DecayMode::BetaMinus),
        gs(55, 137, "Cs-137", 30.04 * yr, DecayMode::BetaMinus),
        gs(56, 133, "Ba-133", 10.5379 * yr, DecayMode::ElectronCapture),
        gs(56, 140, "Ba-140", 12.7534 * dy, DecayMode::BetaMinus),
        gs(57, 140, "La-140", 40.289 * hr, DecayMode::BetaMinus),
        gs(58, 144, "Ce-144", 284.886 * dy, DecayMode::BetaMinus),
        gs(61, 147, "Pm-147", 2.6234 * yr, DecayMode::BetaMinus),
        gs(62, 153, "Sm-153", 46.2846 * hr, DecayMode::BetaMinus),
        gs(63, 152, "Eu-152", 13.517 * yr, DecayMode::ElectronCapture),
        gs(63, 154, "Eu-154", 8.592 * yr, DecayMode::BetaMinus),
        gs(63, 155, "Eu-155", 4.742 * yr, DecayMode::BetaMinus),
        gs(71, 177, "Lu-177", 6.6443 * dy, DecayMode::BetaMinus),
        gs(75, 186, "Re-186", 3.7185 * dy, DecayMode::BetaMinus),
        gs(75, 188, "Re-188", 17.005 * hr, DecayMode::BetaMinus),
        gs(77, 192, "Ir-192", 73.820 * dy, DecayMode::BetaMinus),
        gs(79, 198, "Au-198", 2.69464 * dy, DecayMode::BetaMinus),
        gs(81, 204, "Tl-204", 3.783 * yr, DecayMode::BetaMinus),
        // === Th-232 decay chain ===
        gs(90, 232, "Th-232", 14.0e9 * yr, DecayMode::Alpha),
        gs(88, 228, "Ra-228", 5.75 * yr, DecayMode::BetaMinus),
        gs(89, 228, "Ac-228", 6.15 * hr, DecayMode::BetaMinus),
        gs(90, 228, "Th-228", 1.9125 * yr, DecayMode::Alpha),
        gs(88, 224, "Ra-224", 3.6316 * dy, DecayMode::Alpha),
        gs(86, 220, "Rn-220", 55.6 * sc, DecayMode::Alpha),
        gs(84, 216, "Po-216", 144.0e-3, DecayMode::Alpha),
        gs(82, 212, "Pb-212", 10.627 * hr, DecayMode::BetaMinus),
        gs(83, 212, "Bi-212", 60.55 * mn, DecayMode::BetaMinus), // β⁻ 64.06%, α 35.94%
        gs(81, 208, "Tl-208", 3.053 * mn, DecayMode::BetaMinus),
        gs(84, 212, "Po-212", 294.4e-9, DecayMode::Alpha),
        // === U-235 (actinium) decay chain ===
        gs(90, 231, "Th-231", 25.52 * hr, DecayMode::BetaMinus),
        gs(91, 231, "Pa-231", 32.65e3 * yr, DecayMode::Alpha),
        gs(89, 227, "Ac-227", 21.772 * yr, DecayMode::BetaMinus),
        gs(90, 227, "Th-227", 18.693 * dy, DecayMode::Alpha),
        gs(87, 223, "Fr-223", 22.00 * mn, DecayMode::BetaMinus),
        gs(88, 223, "Ra-223", 11.4352 * dy, DecayMode::Alpha),
        gs(86, 219, "Rn-219", 3.96 * sc, DecayMode::Alpha),
        gs(84, 215, "Po-215", 1.781e-3, DecayMode::Alpha),
        gs(82, 211, "Pb-211", 36.1628 * mn, DecayMode::BetaMinus),
        gs(83, 211, "Bi-211", 2.14 * mn, DecayMode::Alpha),
        gs(81, 207, "Tl-207", 4.77 * mn, DecayMode::BetaMinus),
        // === U-238 decay chain ===
        gs(92, 238, "U-238", 4.463e9 * yr, DecayMode::Alpha),
        gs(90, 234, "Th-234", 24.107 * dy, DecayMode::BetaMinus),
        iso(91, 234, "Pa-234m", 1.159 * mn, DecayMode::BetaMinus, 79.0),
        gs(91, 234, "Pa-234", 6.70 * hr, DecayMode::BetaMinus),
        gs(92, 234, "U-234", 245.5e3 * yr, DecayMode::Alpha),
        gs(92, 235, "U-235", 704.0e6 * yr, DecayMode::Alpha),
        gs(90, 230, "Th-230", 75.4e3 * yr, DecayMode::Alpha),
        gs(88, 226, "Ra-226", 1.600e3 * yr, DecayMode::Alpha),
        gs(86, 222, "Rn-222", 3.8215 * dy, DecayMode::Alpha),
        gs(84, 218, "Po-218", 3.097 * mn, DecayMode::Alpha),
        gs(82, 214, "Pb-214", 27.06 * mn, DecayMode::BetaMinus),
        gs(83, 214, "Bi-214", 19.9 * mn, DecayMode::BetaMinus),
        gs(84, 214, "Po-214", 163.47e-6, DecayMode::Alpha),
        gs(82, 210, "Pb-210", 22.20 * yr, DecayMode::BetaMinus),
        gs(83, 210, "Bi-210", 5.012 * dy, DecayMode::BetaMinus),
        gs(84, 210, "Po-210", 138.376 * dy, DecayMode::Alpha),
        // === Actinides ===
        gs(89, 225, "Ac-225", 9.9190 * dy, DecayMode::Alpha),
        gs(93, 237, "Np-237", 2.144e6 * yr, DecayMode::Alpha),
        gs(94, 238, "Pu-238", 87.7 * yr, DecayMode::Alpha),
        gs(94, 239, "Pu-239", 24.11e3 * yr, DecayMode::Alpha),
        gs(94, 240, "Pu-240", 6.561e3 * yr, DecayMode::Alpha),
        gs(94, 241, "Pu-241", 14.329 * yr, DecayMode::BetaMinus),
        gs(94, 242, "Pu-242", 375.0e3 * yr, DecayMode::Alpha),
        gs(95, 241, "Am-241", 432.6 * yr, DecayMode::Alpha),
        gs(96, 244, "Cm-244", 18.11 * yr, DecayMode::Alpha),
        gs(98, 252, "Cf-252", 2.645 * yr, DecayMode::Alpha),
        // === Isomers ===
        iso(
            43,
            99,
            "Tc-99m",
            6.0066 * hr,
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
        // Ta-180m: no decay has ever been observed; NUBASE2020 gives only the
        // lower limit T½ > 4.5e16 y. The half-life is therefore infinite here and
        // the listed mode is the expected (unobserved) channel.
        iso(
            73,
            180,
            "Ta-180m",
            f64::INFINITY,
            DecayMode::IsomericTransition,
            75.3
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
/// (Bateman solution).
///
/// For a chain N₁ → N₂ → … → N_n with decay constants λ₁ … λ_n and only
/// N₁(0) = `initial_atoms` populated, the populations are the first row of the
/// matrix exponential exp(B), where B is upper bidiagonal with
/// B_ii = −λ_i t and B_i,i+1 = λ_i t:
///
/// N_k(t) = N₁(0) · (Π_{j<k} λ_j t) · exp[−λ₁t, …, −λ_k t],
///
/// the classical Bateman sum written as a divided difference of the
/// exponential. It is evaluated by scaling and squaring with the diagonal and
/// first superdiagonal recomputed exactly after every squaring (Al-Mohy &
/// Higham, SIAM J. Matrix Anal. Appl. 31, 970 (2009); McCurdy, Ng & Parlett,
/// Math. Comp. 43, 501 (1984)). Every intermediate matrix is non-negative, so
/// no cancellation occurs: the result is accurate for equal or nearly equal
/// decay constants, for chains whose λ span many orders of magnitude, and for
/// trace daughters at short times, and it does not depend on the time unit.
///
/// Parameters:
/// - `decay_constants`: λ values in s⁻¹ for each species (0 = stable;
///   negative values are treated as 0)
/// - `initial_atoms`: N₁(0) (only the first species has nonzero initial population)
/// - `time_seconds`: time at which to evaluate (non-positive or non-finite
///   times return the initial state)
///
/// Returns a Vec of atom counts for each species in the chain.
#[must_use]
pub fn bateman_chain(decay_constants: &[f64], initial_atoms: f64, time_seconds: f64) -> Vec<f64> {
    let n = decay_constants.len();
    if n == 0 {
        return Vec::new();
    }
    let mut populations = alloc::vec![0.0; n];
    if time_seconds.is_nan() || time_seconds <= 0.0 || time_seconds.is_infinite() {
        populations[0] = initial_atoms;
        return populations;
    }

    // Scaled nodes x_i = λ_i t ≥ 0.
    let x: Vec<f64> = decay_constants
        .iter()
        .map(|&l| {
            if l > 0.0 && l.is_finite() {
                l * time_seconds
            } else {
                0.0
            }
        })
        .collect();

    let row = exp_chain_first_row(&x);
    for (p, v) in populations.iter_mut().zip(row) {
        *p = initial_atoms * v;
    }
    populations
}

/// Exact superdiagonal entry of exp(h·B) for the 2×2 block
/// [[−a, a], [0, −b]] (with a = λ_i t, b = λ_{i+1} t):
/// h·a·(e^{−ha} − e^{−hb}) / (hb − ha), evaluated without cancellation.
fn chain_superdiag(h: f64, a: f64, b: f64) -> f64 {
    let ha = h * a;
    let hb = h * b;
    let d = 0.5 * (hb - ha);
    if libm::fabs(d) > 0.5 {
        ha * (libm::exp(-ha) - libm::exp(-hb)) / (hb - ha)
    } else {
        // (e^{−ha} − e^{−hb})/(hb − ha) = e^{−(ha+hb)/2} · sinh(d)/d
        let sinhc = if libm::fabs(d) < 1e-8 {
            1.0 + d * d / 6.0
        } else {
            libm::sinh(d) / d
        };
        ha * libm::exp(-0.5 * (ha + hb)) * sinhc
    }
}

/// First row of exp(B) for the chain matrix B (B_ii = −x_i, B_i,i+1 = x_i).
#[allow(clippy::many_single_char_names)]
fn exp_chain_first_row(x: &[f64]) -> Vec<f64> {
    let n = x.len();
    let idx = |i: usize, j: usize| i * n + j;

    // ∞-norm of B is max_i (x_i + x_i) = 2 max x_i. Scale to norm ≤ 1/2.
    let max_x = x.iter().fold(0.0_f64, |m, &v| if v > m { v } else { m });
    let norm = 2.0 * max_x;
    let mut squarings: i32 = 0;
    while libm::ldexp(norm, -squarings) > 0.5 {
        squarings += 1;
    }
    let h = libm::ldexp(1.0, -squarings);

    // Taylor series of exp(h·B). All terms of B^k are exactly representable
    // as products; the norm bound keeps truncation below 1e-19.
    let mut a = alloc::vec![0.0; n * n];
    for i in 0..n {
        a[idx(i, i)] = -h * x[i];
        if i + 1 < n {
            a[idx(i, i + 1)] = h * x[i];
        }
    }
    let mut e = alloc::vec![0.0; n * n];
    for i in 0..n {
        e[idx(i, i)] = 1.0;
    }
    let mut term = e.clone();
    for k in 1..=30 {
        let mut next = alloc::vec![0.0; n * n];
        for i in 0..n {
            for j in i..n {
                let mut acc = 0.0;
                for m in i..=j {
                    acc += term[idx(i, m)] * a[idx(m, j)];
                }
                next[idx(i, j)] = acc / k as f64;
            }
        }
        let mut max_term = 0.0_f64;
        for v in &next {
            max_term = max_term.max(libm::fabs(*v));
        }
        for (ev, tv) in e.iter_mut().zip(next.iter()) {
            *ev += *tv;
        }
        term = next;
        if max_term < 1e-20 {
            break;
        }
    }
    fix_diagonals(&mut e, x, h, n);

    // Square back up, recomputing the diagonal and superdiagonal exactly.
    let mut scale = h;
    for _ in 0..squarings {
        let mut sq = alloc::vec![0.0; n * n];
        for i in 0..n {
            for j in i..n {
                let mut acc = 0.0;
                for m in i..=j {
                    acc += e[idx(i, m)] * e[idx(m, j)];
                }
                sq[idx(i, j)] = acc;
            }
        }
        e = sq;
        scale *= 2.0;
        fix_diagonals(&mut e, x, scale, n);
    }

    (0..n).map(|j| e[idx(0, j)].max(0.0)).collect()
}

/// Overwrites the diagonal and first superdiagonal of exp(scale·B) with their
/// exact values.
fn fix_diagonals(e: &mut [f64], x: &[f64], scale: f64, n: usize) {
    for i in 0..n {
        e[i * n + i] = libm::exp(-scale * x[i]);
        if i + 1 < n {
            e[i * n + i + 1] = chain_superdiag(scale, x[i], x[i + 1]);
        }
    }
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
        let years = c14.half_life_seconds / (365.2422 * 24.0 * 3600.0);
        assert!(
            (years - 5700.0).abs() < 1e-9,
            "C-14 half-life={years} years"
        );
    }

    #[test]
    fn u238_half_life() {
        let u238 = known_isotopes()
            .into_iter()
            .find(|i| i.name == "U-238")
            .unwrap();
        let years = u238.half_life_seconds / (365.2422 * 24.0 * 3600.0);
        assert!(
            (years - 4.463e9).abs() < 1.0,
            "U-238 half-life={years} years"
        );
    }

    #[test]
    fn decay_constant_c14() {
        let t_half = 5700.0 * 365.2422 * 24.0 * 3600.0;
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
    fn bateman_equal_decay_constants() {
        // A → B → C with λ_A = λ_B = 0.1: N_B = N₀ λt e^{−λt}.
        let pops = bateman_chain(&[0.1, 0.1, 0.0], 1.0, 10.0);
        let e = libm::exp(-1.0);
        assert!((pops[0] - e).abs() < 1e-15);
        assert!((pops[1] - e).abs() < 1e-14, "N_B={}", pops[1]);
        assert!((pops[2] - (1.0 - 2.0 * e)).abs() < 1e-14, "N_C={}", pops[2]);
    }

    #[test]
    fn bateman_three_equal_constants() {
        // N_3 = N₀ (λt)²/2 e^{−λt} for three equal λ.
        let pops = bateman_chain(&[0.5, 0.5, 0.5], 1.0, 2.0);
        let expected = 0.5 * libm::exp(-1.0);
        assert!((pops[2] - expected).abs() < 1e-14, "N_3={}", pops[2]);
    }

    #[test]
    fn bateman_independent_of_time_unit() {
        let lam: Vec<f64> = (0..9).map(|i| 1e-13 * (1.0 + f64::from(i))).collect();
        let yr = 3.155_692_6e7;
        let lam_y: Vec<f64> = lam.iter().map(|l| l * yr).collect();
        let a = bateman_chain(&lam, 1.0, 1e13);
        let b = bateman_chain(&lam_y, 1.0, 1e13 / yr);
        // 80-digit reference for the last member: 9.377884e-3.
        assert!((a[8] - 9.377_884e-3).abs() < 1e-9, "last={}", a[8]);
        for (x, y) in a.iter().zip(b.iter()) {
            assert!((x - y).abs() <= 1e-12 * x.abs().max(1e-300), "{x} vs {y}");
        }
    }

    #[test]
    fn bateman_u238_series_trace_daughters() {
        // U-238 series from pure U-238 after 1 year; references computed with
        // the Bateman sum at 80 significant digits.
        let y = 365.2422 * 86400.0;
        let d = 86400.0;
        let m = 60.0;
        let hl = [
            4.468e9 * y,
            24.10 * d,
            1.159 * m,
            2.455e5 * y,
            7.54e4 * y,
            1600.0 * y,
            3.8235 * d,
            3.098 * m,
            27.06 * m,
            19.9 * m,
            164.3e-6,
            22.2 * y,
            5.012 * d,
            138.376 * d,
        ];
        let mut lam: Vec<f64> = hl.iter().map(|h| LN2 / h).collect();
        lam.push(0.0);
        let pops = bateman_chain(&lam, 1.0, y);
        let reference = [
            (5, 5.124_048e-22),  // Ra-226
            (6, 3.193_233e-27),  // Rn-222
            (11, 4.804_590e-26), // Pb-210
            (13, 1.936_077e-28), // Po-210
            (14, 5.539_648e-29), // Pb-206
        ];
        for (i, r) in reference {
            assert!(
                (pops[i] - r).abs() / r < 1e-5,
                "species {i}: {} vs {r}",
                pops[i]
            );
        }
    }

    #[test]
    fn decay_chain_stops_after_isomeric_transition() {
        for (z, a) in [(73, 180), (72, 178), (95, 242)] {
            let chain = decay_chain(&Nucleus::new(z, a).unwrap(), 30);
            let it_steps = chain
                .iter()
                .filter(|(_, m)| *m == DecayMode::IsomericTransition)
                .count();
            assert!(it_steps <= 1, "Z={z} A={a}: {it_steps} IT steps");
        }
    }

    #[test]
    fn bi212_follows_beta_branch() {
        let chain = decay_chain(&Nucleus::new(90, 232).unwrap(), 20);
        let bi = chain
            .iter()
            .find(|(n, _)| n.z() == 83 && n.a() == 212)
            .unwrap();
        assert_eq!(bi.1, DecayMode::BetaMinus);
        assert!(chain.iter().any(|(n, _)| n.z() == 84 && n.a() == 212));
    }

    #[test]
    fn bateman_empty_chain() {
        let pops = bateman_chain(&[], 1000.0, 1.0);
        assert!(pops.is_empty());
    }
}
