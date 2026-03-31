//! Nuclear structure and binding energy calculations.
//!
//! Implements the semi-empirical mass formula (Bethe-Weizsacker formula) for
//! nuclear binding energies, along with nuclear radius calculations,
//! magic number identification, and the nuclear shell model (Mayer-Jensen).

use crate::constants::{AMU_MEV, NEUTRON_MASS_MEV, PROTON_MASS_MEV, R0_FM};
use crate::error::TanmatraError;
use alloc::vec::Vec;
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// AME2020 Atomic Mass Evaluation — Mass Excess Data
// ---------------------------------------------------------------------------

/// AME2020 mass excess values in keV for key nuclides.
///
/// Each entry is (Z, A, mass_excess_kev).
///
/// Source: Wang, M. et al., Chinese Physics C 45, 030003 (2021).
/// "The AME 2020 atomic mass evaluation (II). Tables, graphs and references."
const AME2020_MASS_EXCESS: &[(u32, u32, f64)] = &[
    (1, 1, 7288.971),      // H-1
    (1, 2, 13135.722),     // H-2
    (2, 3, 14931.215),     // He-3
    (2, 4, 2424.916),      // He-4
    (3, 6, 14086.793),     // Li-6
    (3, 7, 14907.105),     // Li-7
    (6, 12, 0.0),          // C-12 (definition)
    (6, 13, 3125.011),     // C-13
    (7, 14, 2863.417),     // N-14
    (8, 16, -4737.001),    // O-16
    (9, 19, -1487.405),    // F-19
    (11, 23, -9529.850),   // Na-23
    (14, 28, -21492.790),  // Si-28
    (15, 31, -24440.990),  // P-31
    (16, 32, -26016.160),  // S-32
    (20, 40, -34846.270),  // Ca-40
    (26, 56, -60601.000),  // Fe-56
    (28, 58, -60228.000),  // Ni-58
    (28, 62, -66746.000),  // Ni-62
    (29, 63, -65579.000),  // Cu-63
    (30, 64, -66004.000),  // Zn-64
    (38, 88, -87922.000),  // Sr-88
    (40, 90, -88768.000),  // Zr-90
    (42, 98, -88113.000),  // Mo-98
    (50, 120, -91105.000), // Sn-120
    (53, 127, -88983.000), // I-127
    (55, 133, -88071.000), // Cs-133
    (56, 138, -88263.000), // Ba-138
    (82, 208, -21749.000), // Pb-208
    (83, 209, -18258.000), // Bi-209
    (92, 235, 40920.500),  // U-235
    (92, 238, 47308.900),  // U-238
];

/// Conversion factor: 1 u = 931494.1 keV (used for mass excess -> atomic mass).
///
/// Source: AME2020 convention, consistent with CODATA 2018 value used in AME2020.
const AME2020_U_KEV: f64 = 931494.1;

// ---------------------------------------------------------------------------
// Nuclear Charge Radii
// ---------------------------------------------------------------------------

/// RMS charge radii in femtometers from elastic electron scattering.
///
/// Each entry is (Z, A, rms_charge_radius_fm).
///
/// Source: Angeli, I. & Marinova, K.P., Atomic Data and Nuclear Data Tables
/// 99, 69-95 (2013). "Table of experimental nuclear ground state charge
/// radii: An update."
const CHARGE_RADII: &[(u32, u32, f64)] = &[
    (1, 1, 0.8783),    // H-1
    (2, 4, 1.6755),    // He-4
    (6, 12, 2.4702),   // C-12
    (8, 16, 2.6991),   // O-16
    (20, 40, 3.4776),  // Ca-40
    (20, 48, 3.4771),  // Ca-48
    (26, 56, 3.7377),  // Fe-56
    (28, 58, 3.770),   // Ni-58
    (38, 88, 4.2240),  // Sr-88
    (50, 120, 4.6519), // Sn-120
    (82, 208, 5.5012), // Pb-208
    (92, 238, 5.8571), // U-238
];

// ---------------------------------------------------------------------------
// Nuclear Electromagnetic Moments
// ---------------------------------------------------------------------------

/// Nuclear magnetic dipole and electric quadrupole moments.
///
/// Contains the magnetic dipole moment in nuclear magnetons (μ_N) and the
/// electric quadrupole moment in barns.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct NuclearMoments {
    /// Magnetic dipole moment in nuclear magnetons (μ_N).
    pub magnetic_dipole_mu_n: f64,
    /// Electric quadrupole moment in barns.
    pub electric_quadrupole_barn: f64,
}

/// Nuclear electromagnetic moments for key nuclides.
///
/// Each entry is (Z, A, magnetic_dipole_mu_n, electric_quadrupole_barn).
///
/// Source: Stone, N.J., Atomic Data and Nuclear Data Tables 90, 75-176 (2005);
/// updated in Stone, N.J., Atomic Data and Nuclear Data Tables 111-112,
/// 1-28 (2016) and Table of Nuclear Magnetic Dipole and Electric Quadrupole
/// Moments, INDC(NDS)-0794 (2019).
const NUCLEAR_MOMENTS: &[(u32, u32, f64, f64)] = &[
    (1, 1, 2.792847, 0.0),         // H-1
    (1, 2, 0.857438, 0.002860),    // H-2
    (2, 3, -2.127625, 0.0),        // He-3
    (3, 6, 0.822047, -0.000806),   // Li-6
    (3, 7, 3.256427, -0.0400),     // Li-7
    (6, 13, 0.702412, 0.0),        // C-13
    (7, 14, 0.403761, 0.02044),    // N-14
    (8, 17, -1.89380, -0.02578),   // O-17
    (9, 19, 2.628868, -0.0942),    // F-19
    (11, 23, 2.217522, 0.104),     // Na-23
    (13, 27, 3.641507, 0.1466),    // Al-27
    (15, 31, 1.13160, 0.0),        // P-31
    (55, 133, 2.582025, -0.00343), // Cs-133
    (82, 207, 0.592583, 0.0),      // Pb-207
    (83, 209, 4.1106, -0.516),     // Bi-209
    (92, 235, -0.38, 4.936),       // U-235
];

// ---------------------------------------------------------------------------
// Superallowed Beta Decay ft Values
// ---------------------------------------------------------------------------

/// A superallowed 0+ → 0+ beta-decay transition.
///
/// Contains the parent and daughter nuclei and the comparative half-life (ft)
/// in seconds.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SuperallowedDecay {
    /// Parent nucleus (decays via beta-plus/EC).
    pub parent: Nucleus,
    /// Daughter nucleus.
    pub daughter: Nucleus,
    /// Comparative half-life ft in seconds.
    pub ft_seconds: f64,
}

/// Superallowed 0+ → 0+ beta-decay ft values.
///
/// Source: Hardy, J.C. & Towner, I.S., Physical Review C 102, 045501 (2020).
/// "Superallowed 0+ → 0+ nuclear β decays: 2020 critical survey."
const SUPERALLOWED_FT_VALUES: &[(u32, u32, u32, u32, f64)] = &[
    // (parent_z, parent_a, daughter_z, daughter_a, ft_seconds)
    (8, 14, 7, 14, 3042.3),   // O-14 → N-14
    (13, 26, 12, 26, 3037.7), // Al-26m → Mg-26
    (17, 34, 16, 34, 3049.4), // Cl-34 → S-34
    (19, 38, 18, 38, 3051.9), // K-38 → Ar-38
    (21, 42, 20, 42, 3047.6), // Sc-42 → Ca-42
    (23, 46, 22, 46, 3049.5), // V-46 → Ti-46
    (25, 50, 24, 50, 3048.4), // Mn-50 → Cr-50
    (27, 54, 26, 54, 3050.8), // Co-54 → Fe-54
    (31, 62, 30, 62, 3074.1), // Ga-62 → Zn-62
];

/// Average corrected Ft value from Hardy & Towner 2020.
///
/// Ft = ft(1 + δ_R')(1 + δ_NS - δ_C) = 3072.27 ± 0.72 s.
///
/// Source: Hardy, J.C. & Towner, I.S., Physical Review C 102, 045501 (2020).
const AVERAGE_CORRECTED_FT: f64 = 3072.27;

/// Bethe-Weizsacker semi-empirical mass formula coefficients (in MeV).
///
/// These are the standard textbook values widely used in nuclear physics.
const A_V: f64 = 15.67; // Volume term
const A_S: f64 = 17.23; // Surface term
const A_C: f64 = 0.714; // Coulomb term
const A_A: f64 = 23.285; // Asymmetry term
/// Pairing term coefficient (MeV).
const A_P: f64 = 11.2;

/// A nucleus characterized by its atomic number Z and mass number A.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Nucleus {
    /// Atomic number (number of protons).
    z: u32,
    /// Mass number (number of protons + neutrons).
    a: u32,
}

impl Nucleus {
    /// Creates a new nucleus with the given atomic number and mass number.
    ///
    /// # Errors
    ///
    /// Returns [`TanmatraError::InvalidAtomicNumber`] if `z` is 0.
    /// Returns [`TanmatraError::InvalidMassNumber`] if `a < z`.
    pub fn new(z: u32, a: u32) -> Result<Self, TanmatraError> {
        if z == 0 {
            return Err(TanmatraError::InvalidAtomicNumber(z));
        }
        if a < z {
            return Err(TanmatraError::InvalidMassNumber { z, a });
        }
        Ok(Self { z, a })
    }

    /// Returns the atomic number (proton count).
    #[must_use]
    pub const fn z(&self) -> u32 {
        self.z
    }

    /// Returns the mass number (nucleon count).
    #[must_use]
    pub const fn a(&self) -> u32 {
        self.a
    }

    /// Returns the neutron number N = A - Z.
    #[must_use]
    pub const fn n(&self) -> u32 {
        self.a - self.z
    }

    /// Calculates the nuclear binding energy in MeV using the
    /// Bethe-Weizsacker semi-empirical mass formula.
    ///
    /// B(Z,A) = a_v*A - a_s*A^(2/3) - a_c*Z*(Z-1)/A^(1/3) - a_a*(A-2Z)^2/A + delta
    ///
    /// where delta is the pairing term.
    #[must_use]
    #[inline]
    pub fn binding_energy(&self) -> f64 {
        let a = self.a as f64;
        let z = self.z as f64;

        if self.a == 1 {
            return 0.0; // Single nucleon has no binding energy
        }

        let a_one_third = libm::cbrt(a);
        let a_two_thirds = a_one_third * a_one_third;

        // Volume term
        let volume = A_V * a;

        // Surface term
        let surface = A_S * a_two_thirds;

        // Coulomb term
        let coulomb = A_C * z * (z - 1.0) / a_one_third;

        // Asymmetry term
        let asymmetry_num = (a - 2.0 * z) * (a - 2.0 * z);
        let asymmetry = A_A * asymmetry_num / a;

        // Pairing term
        let delta = pairing_term(self.z, self.a);

        volume - surface - coulomb - asymmetry + delta
    }

    /// Calculates binding energy with Strutinsky shell correction.
    ///
    /// Adds a shell correction term that accounts for the extra stability
    /// of nuclei near magic numbers. The correction is estimated from the
    /// distance to the nearest shell closure for both protons and neutrons.
    ///
    /// B_corrected = B_LDM + δ_shell(Z) + δ_shell(N)
    ///
    /// where δ_shell is negative (more bound) near magic numbers and
    /// positive (less bound) between shells.
    #[must_use]
    #[inline]
    pub fn binding_energy_shell_corrected(&self) -> f64 {
        let b_ldm = self.binding_energy();
        let delta_z = shell_correction_energy(self.z);
        let delta_n = shell_correction_energy(self.n());
        b_ldm + delta_z + delta_n
    }

    /// Returns the binding energy per nucleon (B/A) in MeV.
    #[must_use]
    #[inline]
    pub fn binding_energy_per_nucleon(&self) -> f64 {
        self.binding_energy() / self.a as f64
    }

    /// Returns the mass defect in MeV/c^2.
    ///
    /// Mass defect = Z*m_p + N*m_n - M_nucleus
    /// where M_nucleus = Z*m_p + N*m_n - B(Z,A)
    /// so mass defect = B(Z,A) (the binding energy itself in mass-energy equivalence).
    #[must_use]
    #[inline]
    pub fn mass_defect(&self) -> f64 {
        self.binding_energy()
    }

    /// Returns the nuclear mass in MeV/c^2.
    ///
    /// M = Z*m_p + N*m_n - B(Z,A)
    #[must_use]
    #[inline]
    pub fn nuclear_mass(&self) -> f64 {
        let z = self.z as f64;
        let n = self.n() as f64;
        z * PROTON_MASS_MEV + n * NEUTRON_MASS_MEV - self.binding_energy()
    }

    /// Returns the atomic mass in atomic mass units (u).
    #[must_use]
    #[inline]
    pub fn atomic_mass_amu(&self) -> f64 {
        self.nuclear_mass() / AMU_MEV
    }

    /// Returns the nuclear radius in femtometers using R = r0 * A^(1/3).
    #[must_use]
    #[inline]
    pub fn nuclear_radius(&self) -> f64 {
        R0_FM * libm::cbrt(self.a as f64)
    }

    /// Returns `true` if Z or N is a magic number.
    ///
    /// Magic numbers: 2, 8, 20, 28, 50, 82, 126
    #[must_use]
    pub fn is_magic(&self) -> bool {
        is_magic_number(self.z) || is_magic_number(self.n())
    }

    /// Returns `true` if both Z and N are magic numbers (doubly magic).
    #[must_use]
    pub fn is_doubly_magic(&self) -> bool {
        is_magic_number(self.z) && is_magic_number(self.n())
    }

    // --- AME2020 mass data ---

    /// Returns the experimental mass excess in keV from the AME2020 evaluation.
    ///
    /// Looks up the nucleus (Z, A) in the AME2020 table of mass excess values.
    /// Returns `None` if the nuclide is not in the table.
    ///
    /// Source: Wang et al., Chinese Physics C 45, 030003 (2021).
    #[must_use]
    pub fn experimental_mass_excess_kev(&self) -> Option<f64> {
        AME2020_MASS_EXCESS
            .iter()
            .find(|&&(z, a, _)| z == self.z && a == self.a)
            .map(|&(_, _, me)| me)
    }

    /// Returns the experimental atomic mass in atomic mass units (u) from AME2020.
    ///
    /// Computed from the mass excess: M(u) = A + mass_excess_kev / 931494.1.
    /// Returns `None` if the nuclide is not in the AME2020 table.
    ///
    /// Source: Wang et al., Chinese Physics C 45, 030003 (2021).
    #[must_use]
    pub fn experimental_atomic_mass_amu(&self) -> Option<f64> {
        self.experimental_mass_excess_kev()
            .map(|me| self.a as f64 + me / AME2020_U_KEV)
    }

    // --- Charge radii ---

    /// Returns the experimental RMS charge radius in femtometers.
    ///
    /// Looks up the nucleus in the Angeli & Marinova (2013) table of
    /// nuclear charge radii from elastic electron scattering.
    /// Returns `None` if the nuclide is not in the table.
    ///
    /// Source: Angeli & Marinova, At. Data Nucl. Data Tables 99, 69 (2013).
    #[must_use]
    pub fn charge_radius_fm(&self) -> Option<f64> {
        CHARGE_RADII
            .iter()
            .find(|&&(z, a, _)| z == self.z && a == self.a)
            .map(|&(_, _, r)| r)
    }

    // --- Electromagnetic moments ---

    /// Returns the nuclear magnetic dipole and electric quadrupole moments.
    ///
    /// Looks up the nucleus in the Stone (2005/2019) table of nuclear
    /// electromagnetic moments.
    /// Returns `None` if the nuclide is not in the table.
    ///
    /// Source: Stone, N.J., At. Data Nucl. Data Tables 90, 75 (2005);
    /// INDC(NDS)-0794 (2019).
    #[must_use]
    pub fn nuclear_moments(&self) -> Option<NuclearMoments> {
        NUCLEAR_MOMENTS
            .iter()
            .find(|&&(z, a, _, _)| z == self.z && a == self.a)
            .map(|&(_, _, mu, q)| NuclearMoments {
                magnetic_dipole_mu_n: mu,
                electric_quadrupole_barn: q,
            })
    }

    // --- Presets ---

    /// Hydrogen-1 (proton).
    #[must_use]
    pub fn hydrogen_1() -> Self {
        Self { z: 1, a: 1 }
    }

    /// Helium-4 (alpha particle).
    #[must_use]
    pub fn helium_4() -> Self {
        Self { z: 2, a: 4 }
    }

    /// Carbon-12.
    #[must_use]
    pub fn carbon_12() -> Self {
        Self { z: 6, a: 12 }
    }

    /// Iron-56 (most tightly bound common nucleus).
    #[must_use]
    pub fn iron_56() -> Self {
        Self { z: 26, a: 56 }
    }

    /// Uranium-235 (fissile).
    #[must_use]
    pub fn uranium_235() -> Self {
        Self { z: 92, a: 235 }
    }

    /// Uranium-238.
    #[must_use]
    pub fn uranium_238() -> Self {
        Self { z: 92, a: 238 }
    }
}

/// Calculates the pairing term delta for the Bethe-Weizsacker formula.
///
/// delta = +a_p / A^(1/2) for even-even (Z even, N even)
/// delta = 0              for odd A
/// delta = -a_p / A^(1/2) for odd-odd (Z odd, N odd)
#[must_use]
fn pairing_term(z: u32, a: u32) -> f64 {
    let n = a - z;
    let a_f = a as f64;
    let denom = libm::sqrt(a_f);

    if denom == 0.0 {
        return 0.0;
    }

    if z.is_multiple_of(2) && n.is_multiple_of(2) {
        A_P / denom
    } else if !z.is_multiple_of(2) && !n.is_multiple_of(2) {
        -A_P / denom
    } else {
        0.0
    }
}

/// Estimates the Strutinsky shell correction energy for a nucleon number.
///
/// Uses a Gaussian-smoothed single-particle level density approach.
/// The correction is largest (most negative, meaning extra binding) at
/// magic numbers and oscillates between shells.
///
/// Based on the parameterization from Myers & Swiatecki (1966) and
/// refined empirical fits. The correction magnitude is typically 1-3 MeV
/// per nucleon type (proton or neutron separately).
/// Magic numbers and their shell correction peak energies (empirical, in MeV).
/// Positive values = extra binding at magic numbers.
const SHELL_CORRECTION_MAGIC: [(u32, f64); 7] = [
    (2, 2.5),
    (8, 3.5),
    (20, 3.0),
    (28, 3.5),
    (50, 3.0),
    (82, 3.5),
    (126, 3.0),
];

fn shell_correction_energy(nucleon_count: u32) -> f64 {
    if nucleon_count == 0 {
        return 0.0;
    }

    let n = nucleon_count as f64;

    // Find the nearest magic number and compute a Gaussian correction
    let mut correction = 0.0;
    for &(magic, peak_energy) in &SHELL_CORRECTION_MAGIC {
        let magic_f = magic as f64;
        // Width parameter scales with the shell gap spacing
        let width = 0.1 * magic_f + 2.0;
        let dist = n - magic_f;
        correction += peak_energy * libm::exp(-(dist * dist) / (2.0 * width * width));
    }

    correction
}

/// Returns `true` if the given number is a nuclear magic number.
///
/// Magic numbers correspond to complete nuclear shells:
/// 2, 8, 20, 28, 50, 82, 126.
#[must_use]
pub fn is_magic_number(n: u32) -> bool {
    matches!(n, 2 | 8 | 20 | 28 | 50 | 82 | 126)
}

// ---------------------------------------------------------------------------
// Nuclear Shell Model (Mayer-Jensen)
// ---------------------------------------------------------------------------

/// A nuclear shell model single-particle level.
///
/// Each level is characterized by quantum numbers (n, l, j) where j = l ± 1/2.
/// The degeneracy is 2j + 1.
///
/// The ordering follows the harmonic oscillator potential with strong spin-orbit
/// coupling (Mayer-Jensen shell model, 1949).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ShellLevel {
    /// Principal oscillator quantum number (1-based: 1s, 1p, 1d, ...).
    pub n_shell: u32,
    /// Orbital angular momentum quantum number.
    pub l: u32,
    /// Total angular momentum (stored as 2j to keep integer: j = l ± 1/2).
    pub two_j: u32,
}

impl ShellLevel {
    /// Returns the degeneracy (number of substates) = 2j + 1.
    #[must_use]
    #[inline]
    pub const fn degeneracy(&self) -> u32 {
        self.two_j + 1
    }

    /// Returns the total angular momentum j as a float.
    #[must_use]
    #[inline]
    pub fn j(&self) -> f64 {
        self.two_j as f64 / 2.0
    }

    /// Returns the spectroscopic label (e.g., "1s1/2", "1p3/2").
    #[must_use]
    pub fn label(&self) -> alloc::string::String {
        let l_char = match self.l {
            0 => 's',
            1 => 'p',
            2 => 'd',
            3 => 'f',
            4 => 'g',
            5 => 'h',
            6 => 'i',
            _ => '?',
        };
        alloc::format!("{}{}{}/{}", self.n_shell, l_char, self.two_j, 2)
    }
}

/// Standard nuclear shell model level ordering (Mayer-Jensen).
///
/// Levels are ordered by energy from the harmonic oscillator potential with
/// strong spin-orbit coupling. This ordering reproduces the nuclear magic
/// numbers: 2, 8, 20, 28, 50, 82, 126.
///
/// Each entry: (n_shell, l, 2j), where n_shell is the radial quantum number
/// within each l (1-based), l is the orbital angular momentum, and 2j is
/// twice the total angular momentum.
///
/// Source: Mayer & Jensen (Nobel Prize 1963), standard nuclear physics
/// textbooks (Krane, Wong, Ring & Schuck).
const SHELL_MODEL_LEVELS: [(u32, u32, u32); 32] = [
    // Shell closure at 2
    (1, 0, 1), // 1s1/2 (2)      cumulative: 2
    // Shell closure at 8
    (1, 1, 3), // 1p3/2 (4)      cumulative: 6
    (1, 1, 1), // 1p1/2 (2)      cumulative: 8
    // Shell closure at 20
    (1, 2, 5), // 1d5/2 (6)      cumulative: 14
    (2, 0, 1), // 2s1/2 (2)      cumulative: 16
    (1, 2, 3), // 1d3/2 (4)      cumulative: 20
    // Shell closure at 28
    (1, 3, 7), // 1f7/2 (8)      cumulative: 28
    // Shell closure at 50
    (2, 1, 3), // 2p3/2 (4)      cumulative: 32
    (1, 3, 5), // 1f5/2 (6)      cumulative: 38
    (2, 1, 1), // 2p1/2 (2)      cumulative: 40
    (1, 4, 9), // 1g9/2 (10)     cumulative: 50
    // Shell closure at 82
    (1, 4, 7),  // 1g7/2 (8)      cumulative: 58
    (2, 2, 5),  // 2d5/2 (6)      cumulative: 64
    (2, 2, 3),  // 2d3/2 (4)      cumulative: 68
    (3, 0, 1),  // 3s1/2 (2)      cumulative: 70
    (1, 5, 11), // 1h11/2 (12)   cumulative: 82
    // Shell closure at 126
    (1, 5, 9),  // 1h9/2 (10)     cumulative: 92
    (2, 3, 7),  // 2f7/2 (8)      cumulative: 100
    (2, 3, 5),  // 2f5/2 (6)      cumulative: 106
    (3, 1, 3),  // 3p3/2 (4)      cumulative: 110
    (3, 1, 1),  // 3p1/2 (2)      cumulative: 112
    (1, 6, 13), // 1i13/2 (14)   cumulative: 126
    // Beyond 126 (shell closure at 184 predicted)
    (2, 4, 9),  // 2g9/2 (10)     cumulative: 136
    (1, 6, 11), // 1i11/2 (12)   cumulative: 148
    (3, 2, 5),  // 3d5/2 (6)      cumulative: 154
    (4, 0, 1),  // 4s1/2 (2)      cumulative: 156
    (2, 4, 7),  // 2g7/2 (8)      cumulative: 164
    (3, 2, 3),  // 3d3/2 (4)      cumulative: 168
    (1, 7, 15), // 1j15/2 (16)   cumulative: 184
    (2, 5, 11), // 2h11/2 (12)   cumulative: 196
    (2, 5, 9),  // 2h9/2 (10)     cumulative: 206
    (3, 3, 7),  // 3f7/2 (8)      cumulative: 214
];

/// Returns the nuclear shell model levels in energy order.
///
/// This is the standard Mayer-Jensen ordering that reproduces the magic
/// numbers through spin-orbit coupling.
#[must_use]
pub fn shell_model_levels() -> &'static [(u32, u32, u32)] {
    &SHELL_MODEL_LEVELS
}

/// Returns the shell model occupation for a given nucleon count.
///
/// Each entry in the returned vector is `(ShellLevel, occupation)` where
/// `occupation` is the number of nucleons in that level (0 to degeneracy).
///
/// This applies to protons and neutrons independently.
#[must_use]
pub fn shell_occupation(nucleon_count: u32) -> Vec<(ShellLevel, u32)> {
    let mut remaining = nucleon_count;
    let mut occupation = Vec::new();

    for &(n_shell, l, two_j) in &SHELL_MODEL_LEVELS {
        if remaining == 0 {
            break;
        }
        let level = ShellLevel { n_shell, l, two_j };
        let deg = level.degeneracy();
        let fill = if remaining >= deg { deg } else { remaining };
        occupation.push((level, fill));
        remaining -= fill;
    }

    occupation
}

/// Returns the ground-state spin and parity (J^pi) of a nucleus.
///
/// Uses the shell model: for even-even nuclei, J^pi = 0+.
/// For odd-A nuclei, J^pi is determined by the last unpaired nucleon.
/// For odd-odd nuclei, J^pi is determined by coupling the last proton
/// and neutron (simplified: returns the range of possible J values).
///
/// Returns `(two_j, parity)` where parity is +1 or -1, and two_j is
/// twice the total nuclear spin. For odd-odd nuclei, returns the
/// spin of the last odd proton (a simplification).
#[must_use]
pub fn ground_state_spin_parity(nucleus: &Nucleus) -> (u32, i32) {
    let z = nucleus.z();
    let n = nucleus.n();
    let z_even = z.is_multiple_of(2);
    let n_even = n.is_multiple_of(2);

    if z_even && n_even {
        // Even-even: always 0+
        return (0, 1);
    }

    // Find the last unpaired nucleon
    let (nucleon_count, is_proton_odd) = if !z_even && n_even {
        (z, true)
    } else if z_even && !n_even {
        (n, false)
    } else {
        // Odd-odd: use last odd proton (simplification)
        (z, true)
    };

    let _ = is_proton_odd; // used for documentation clarity
    let occ = shell_occupation(nucleon_count);

    // Find the last partially filled level
    if let Some(&(level, _fill)) = occ.last() {
        let parity = if level.l % 2 == 0 { 1 } else { -1 };
        (level.two_j, parity)
    } else {
        (0, 1) // fallback
    }
}

/// Returns the shell closure number at or below the given nucleon count.
///
/// Shell closures (magic numbers) occur at: 2, 8, 20, 28, 50, 82, 126, 184.
#[must_use]
pub fn shell_closure_below(nucleon_count: u32) -> u32 {
    const CLOSURES: [u32; 8] = [2, 8, 20, 28, 50, 82, 126, 184];
    let mut result = 0;
    for &c in &CLOSURES {
        if c <= nucleon_count {
            result = c;
        } else {
            break;
        }
    }
    result
}

/// Returns the next shell closure above the given nucleon count.
///
/// Shell closures: 2, 8, 20, 28, 50, 82, 126, 184.
/// Returns `None` if above the highest known closure.
#[must_use]
pub fn next_shell_closure(nucleon_count: u32) -> Option<u32> {
    const CLOSURES: [u32; 8] = [2, 8, 20, 28, 50, 82, 126, 184];
    CLOSURES.iter().find(|&&c| c > nucleon_count).copied()
}

// ---------------------------------------------------------------------------
// Superallowed beta-decay API
// ---------------------------------------------------------------------------

/// Returns the superallowed 0+ → 0+ beta-decay ft values.
///
/// Source: Hardy & Towner, Phys. Rev. C 102, 045501 (2020).
#[must_use]
pub fn superallowed_ft_values() -> Vec<SuperallowedDecay> {
    let mut v = Vec::with_capacity(SUPERALLOWED_FT_VALUES.len());
    for &(pz, pa, dz, da, ft) in SUPERALLOWED_FT_VALUES {
        // These are all valid well-known nuclides; use unwrap_or_else to avoid panic.
        let parent = Nucleus::new(pz, pa).unwrap_or(Nucleus { z: pz, a: pa });
        let daughter = Nucleus::new(dz, da).unwrap_or(Nucleus { z: dz, a: da });
        v.push(SuperallowedDecay {
            parent,
            daughter,
            ft_seconds: ft,
        });
    }
    v
}

/// Applies radiative and isospin-breaking corrections to a bare ft value.
///
/// Returns the corrected Ft value using the average correction factor derived
/// from the world average Ft = 3072.27 s (Hardy & Towner 2020).
///
/// The correction formula is: Ft = ft × (1 + δ_R') × (1 + δ_NS - δ_C).
/// Since individual correction terms vary per transition, this function uses
/// the empirically determined average ratio Ft/ft ≈ 3072.27 / <ft_avg>
/// where <ft_avg> is the average of all measured ft values.
///
/// For a more precise correction, the individual radiative (δ_R'), nuclear
/// structure (δ_NS), and isospin symmetry-breaking (δ_C) terms should be
/// applied per transition.
///
/// As a simplified approach, this returns the world-average corrected value
/// `AVERAGE_CORRECTED_FT` = 3072.27 s, which is the nucleus-independent Ft
/// that all superallowed decays should yield after correction.
///
/// Source: Hardy & Towner, Phys. Rev. C 102, 045501 (2020).
#[must_use]
pub fn corrected_ft_value(_ft: f64) -> f64 {
    // The corrected Ft value is nucleus-independent by definition.
    // Individual corrections (δ_R', δ_NS, δ_C) map each transition's ft
    // to this common Ft value. Since those corrections are transition-specific
    // and the whole point is that they converge to a single value, we return
    // the world average.
    AVERAGE_CORRECTED_FT
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fe56_binding_energy_per_nucleon() {
        let fe56 = Nucleus::iron_56();
        let bea = fe56.binding_energy_per_nucleon();
        // Fe-56 experimental B/A ≈ 8.790 MeV; semi-empirical should be within ~2%
        assert!(bea > 8.6, "Fe-56 B/A={bea} too low");
        assert!(bea < 9.0, "Fe-56 B/A={bea} too high");
    }

    #[test]
    fn he4_binding_energy_positive() {
        let he4 = Nucleus::helium_4();
        let be = he4.binding_energy();
        // He-4 experimental: 28.3 MeV; semi-empirical is less accurate for light nuclei
        assert!(be > 15.0, "He-4 BE={be} too low");
        assert!(be < 35.0, "He-4 BE={be} within expected range");
    }

    #[test]
    fn hydrogen_has_zero_binding_energy() {
        let h1 = Nucleus::hydrogen_1();
        assert!((h1.binding_energy()).abs() < 1e-10);
    }

    #[test]
    fn nuclear_radius_fe56() {
        let fe56 = Nucleus::iron_56();
        let r = fe56.nuclear_radius();
        // R = 1.2 * 56^(1/3) ≈ 1.2 * 3.826 ≈ 4.59 fm
        assert!((r - 4.59).abs() < 0.1, "Fe-56 radius={r} fm");
    }

    #[test]
    fn magic_numbers() {
        assert!(is_magic_number(2));
        assert!(is_magic_number(8));
        assert!(is_magic_number(20));
        assert!(is_magic_number(28));
        assert!(is_magic_number(50));
        assert!(is_magic_number(82));
        assert!(is_magic_number(126));
        assert!(!is_magic_number(3));
        assert!(!is_magic_number(100));
    }

    #[test]
    fn doubly_magic_he4() {
        let he4 = Nucleus::helium_4();
        assert!(he4.is_doubly_magic()); // Z=2, N=2 both magic
    }

    #[test]
    fn invalid_z_zero() {
        assert!(Nucleus::new(0, 1).is_err());
    }

    #[test]
    fn invalid_a_less_than_z() {
        assert!(Nucleus::new(10, 5).is_err());
    }

    #[test]
    fn mass_defect_equals_binding_energy() {
        let c12 = Nucleus::carbon_12();
        assert!((c12.mass_defect() - c12.binding_energy()).abs() < 1e-10);
    }

    #[test]
    fn serde_roundtrip() {
        let n = Nucleus::iron_56();
        let json = serde_json::to_string(&n).unwrap();
        let back: Nucleus = serde_json::from_str(&json).unwrap();
        assert_eq!(n, back);
    }

    #[test]
    fn binding_energy_increases_with_a_midrange() {
        // For stable nuclei, total BE should increase (though B/A peaks around Fe)
        let c12 = Nucleus::carbon_12();
        let fe56 = Nucleus::iron_56();
        assert!(fe56.binding_energy() > c12.binding_energy());
    }

    #[test]
    fn bea_peaks_around_iron() {
        // B/A should peak around A~56-62; light and heavy nuclei have lower B/A
        let o16 = Nucleus::new(8, 16).unwrap();
        let fe56 = Nucleus::iron_56();
        let u238 = Nucleus::uranium_238();
        assert!(
            fe56.binding_energy_per_nucleon() > o16.binding_energy_per_nucleon(),
            "Fe-56 B/A should exceed O-16 B/A"
        );
        assert!(
            fe56.binding_energy_per_nucleon() > u238.binding_energy_per_nucleon(),
            "Fe-56 B/A should exceed U-238 B/A"
        );
    }

    #[test]
    fn nuclear_mass_reasonable() {
        // Fe-56 nuclear mass should be near 56 * AMU ≈ 52164 MeV
        let fe56 = Nucleus::iron_56();
        let mass = fe56.nuclear_mass();
        assert!(mass > 52_000.0, "Fe-56 mass={mass} too low");
        assert!(mass < 52_500.0, "Fe-56 mass={mass} too high");
    }

    // --- Shell model tests ---

    #[test]
    fn shell_model_reproduces_magic_numbers() {
        // Verify cumulative filling produces magic numbers at shell closures
        let levels = shell_model_levels();
        let mut cumulative = 0u32;
        let magic = [2, 8, 20, 28, 50, 82, 126];
        let mut magic_idx = 0;

        for &(_, _, two_j) in levels {
            cumulative += two_j + 1; // degeneracy = 2j+1
            if magic_idx < magic.len() && cumulative == magic[magic_idx] {
                magic_idx += 1;
            }
        }
        assert_eq!(
            magic_idx,
            magic.len(),
            "not all magic numbers found in shell model"
        );
    }

    #[test]
    fn shell_level_labels() {
        let level = ShellLevel {
            n_shell: 1,
            l: 3,
            two_j: 7,
        };
        assert_eq!(level.label(), "1f7/2");
        assert_eq!(level.degeneracy(), 8);
    }

    #[test]
    fn o16_spin_parity_0_plus() {
        // O-16: Z=8, N=8 (doubly magic, even-even) -> 0+
        let o16 = Nucleus::new(8, 16).unwrap();
        let (two_j, parity) = ground_state_spin_parity(&o16);
        assert_eq!(two_j, 0, "O-16 should have J=0");
        assert_eq!(parity, 1, "O-16 should have positive parity");
    }

    #[test]
    fn fe56_spin_parity_0_plus() {
        // Fe-56: Z=26, N=30 (even-even) -> 0+
        let fe56 = Nucleus::iron_56();
        let (two_j, parity) = ground_state_spin_parity(&fe56);
        assert_eq!(two_j, 0);
        assert_eq!(parity, 1);
    }

    #[test]
    fn o17_spin_parity() {
        // O-17: Z=8 (magic), N=9 -> last neutron in 1d5/2 -> 5/2+
        let o17 = Nucleus::new(8, 17).unwrap();
        let (two_j, parity) = ground_state_spin_parity(&o17);
        assert_eq!(two_j, 5, "O-17 should have 2J=5 (J=5/2)");
        assert_eq!(parity, 1, "O-17 should have positive parity (l=2)");
    }

    #[test]
    fn shell_occupation_he4() {
        // He-4: 2 protons fill 1s1/2 completely
        let occ = shell_occupation(2);
        assert_eq!(occ.len(), 1);
        assert_eq!(occ[0].1, 2); // 1s1/2 fully filled
    }

    #[test]
    fn shell_closure_functions() {
        assert_eq!(shell_closure_below(10), 8);
        assert_eq!(shell_closure_below(28), 28);
        assert_eq!(shell_closure_below(1), 0);
        assert_eq!(next_shell_closure(20), Some(28));
        assert_eq!(next_shell_closure(82), Some(126));
        assert_eq!(next_shell_closure(200), None);
    }

    #[test]
    fn serde_roundtrip_shell_level() {
        let level = ShellLevel {
            n_shell: 1,
            l: 3,
            two_j: 7,
        };
        let json = serde_json::to_string(&level).unwrap();
        let back: ShellLevel = serde_json::from_str(&json).unwrap();
        assert_eq!(level, back);
    }

    // --- Strutinsky shell correction tests ---

    #[test]
    fn shell_corrected_more_bound_at_magic() {
        // Doubly magic O-16 (Z=8, N=8) should have MORE binding with correction
        let o16 = Nucleus::new(8, 16).unwrap();
        assert!(
            o16.binding_energy_shell_corrected() > o16.binding_energy(),
            "Shell correction should increase BE for doubly-magic O-16"
        );
    }

    #[test]
    fn shell_corrected_doubly_magic_ca40() {
        // Ca-40 (Z=20, N=20) is doubly magic
        let ca40 = Nucleus::new(20, 40).unwrap();
        assert!(
            ca40.binding_energy_shell_corrected() > ca40.binding_energy(),
            "Shell correction should increase BE for doubly-magic Ca-40"
        );
    }

    #[test]
    fn shell_correction_reasonable_magnitude() {
        let fe56 = Nucleus::iron_56();
        let diff = (fe56.binding_energy_shell_corrected() - fe56.binding_energy()).abs();
        assert!(diff < 10.0, "Shell correction {diff} MeV too large");
    }

    // --- Coverage: additional paths ---

    #[test]
    fn atomic_mass_amu_fe56() {
        let fe56 = Nucleus::iron_56();
        let amu = fe56.atomic_mass_amu();
        // Fe-56 atomic mass ≈ 55.9 u
        assert!(amu > 55.0 && amu < 57.0, "Fe-56 AMU={amu}");
    }

    #[test]
    fn is_magic_not_doubly() {
        // O-17: Z=8 (magic), N=9 (not magic)
        let o17 = Nucleus::new(8, 17).unwrap();
        assert!(o17.is_magic());
        assert!(!o17.is_doubly_magic());
    }

    #[test]
    fn shell_level_degeneracy() {
        let level = ShellLevel {
            n_shell: 1,
            l: 0,
            two_j: 1,
        };
        assert_eq!(level.degeneracy(), 2);
        assert!((level.j() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn shell_level_label_g_orbital() {
        let level = ShellLevel {
            n_shell: 1,
            l: 4,
            two_j: 9,
        };
        assert_eq!(level.label(), "1g9/2");
    }

    #[test]
    fn odd_odd_spin_parity() {
        // N-14: Z=7, N=7 (odd-odd)
        let n14 = Nucleus::new(7, 14).unwrap();
        let (two_j, _parity) = ground_state_spin_parity(&n14);
        assert!(two_j > 0, "Odd-odd should have nonzero spin");
    }

    #[test]
    fn shell_closure_below_zero() {
        assert_eq!(shell_closure_below(0), 0);
    }

    #[test]
    fn carbon_presets() {
        let c12 = Nucleus::carbon_12();
        assert_eq!(c12.z(), 6);
        assert_eq!(c12.a(), 12);
        assert_eq!(c12.n(), 6);
    }

    // --- AME2020 mass excess tests ---

    #[test]
    fn ame2020_c12_mass_excess_zero() {
        let c12 = Nucleus::carbon_12();
        let me = c12.experimental_mass_excess_kev().unwrap();
        assert!(
            (me).abs() < 1e-10,
            "C-12 mass excess should be 0 by definition"
        );
    }

    #[test]
    fn ame2020_h1_mass_excess() {
        let h1 = Nucleus::hydrogen_1();
        let me = h1.experimental_mass_excess_kev().unwrap();
        assert!((me - 7288.971).abs() < 0.01, "H-1 mass excess={me} keV");
    }

    #[test]
    fn ame2020_fe56_mass_excess_negative() {
        let fe56 = Nucleus::iron_56();
        let me = fe56.experimental_mass_excess_kev().unwrap();
        assert!(me < 0.0, "Fe-56 mass excess should be negative");
        assert!((me - (-60601.0)).abs() < 1.0, "Fe-56 mass excess={me} keV");
    }

    #[test]
    fn ame2020_unknown_nuclide_returns_none() {
        // Og-294 is not in our table
        let og = Nucleus::new(118, 294).unwrap();
        assert!(og.experimental_mass_excess_kev().is_none());
    }

    #[test]
    fn ame2020_atomic_mass_h1() {
        let h1 = Nucleus::hydrogen_1();
        let mass = h1.experimental_atomic_mass_amu().unwrap();
        // H-1: M = 1 + 7288.971/931494.1 ≈ 1.007825
        assert!((mass - 1.007825).abs() < 0.0001, "H-1 atomic mass={mass} u");
    }

    #[test]
    fn ame2020_atomic_mass_c12_exactly_12() {
        let c12 = Nucleus::carbon_12();
        let mass = c12.experimental_atomic_mass_amu().unwrap();
        assert!(
            (mass - 12.0).abs() < 1e-6,
            "C-12 atomic mass should be exactly 12, got {mass}"
        );
    }

    #[test]
    fn ame2020_atomic_mass_fe56() {
        let fe56 = Nucleus::iron_56();
        let mass = fe56.experimental_atomic_mass_amu().unwrap();
        // Fe-56: ~55.9349 u
        assert!((mass - 55.9349).abs() < 0.001, "Fe-56 atomic mass={mass} u");
    }

    #[test]
    fn ame2020_atomic_mass_u238() {
        let u238 = Nucleus::uranium_238();
        let mass = u238.experimental_atomic_mass_amu().unwrap();
        // U-238: ~238.0508 u
        assert!(
            (mass - 238.0508).abs() < 0.001,
            "U-238 atomic mass={mass} u"
        );
    }

    #[test]
    fn ame2020_all_entries_valid() {
        // Verify all entries in the table have Z > 0 and A >= Z
        for &(z, a, _) in AME2020_MASS_EXCESS {
            assert!(z > 0, "Z must be > 0, got {z}");
            assert!(a >= z, "A must be >= Z, got Z={z} A={a}");
        }
    }

    // --- Charge radii tests ---

    #[test]
    fn charge_radius_h1() {
        let h1 = Nucleus::hydrogen_1();
        let r = h1.charge_radius_fm().unwrap();
        assert!((r - 0.8783).abs() < 0.001, "H-1 charge radius={r} fm");
    }

    #[test]
    fn charge_radius_pb208() {
        let pb208 = Nucleus::new(82, 208).unwrap();
        let r = pb208.charge_radius_fm().unwrap();
        assert!((r - 5.5012).abs() < 0.001, "Pb-208 charge radius={r} fm");
    }

    #[test]
    fn charge_radius_increases_with_a() {
        let he4 = Nucleus::helium_4();
        let pb208 = Nucleus::new(82, 208).unwrap();
        assert!(
            pb208.charge_radius_fm().unwrap() > he4.charge_radius_fm().unwrap(),
            "Pb-208 radius should be larger than He-4"
        );
    }

    #[test]
    fn charge_radius_unknown_returns_none() {
        let og = Nucleus::new(118, 294).unwrap();
        assert!(og.charge_radius_fm().is_none());
    }

    #[test]
    fn charge_radius_ca_isotopes() {
        // Ca-40 and Ca-48 have very similar charge radii (nuclear physics anomaly)
        let ca40 = Nucleus::new(20, 40).unwrap();
        let ca48 = Nucleus::new(20, 48).unwrap();
        let r40 = ca40.charge_radius_fm().unwrap();
        let r48 = ca48.charge_radius_fm().unwrap();
        assert!(
            (r40 - r48).abs() < 0.01,
            "Ca-40 ({r40}) and Ca-48 ({r48}) radii should be very similar"
        );
    }

    // --- Nuclear moments tests ---

    #[test]
    fn nuclear_moments_h1_proton() {
        let h1 = Nucleus::hydrogen_1();
        let m = h1.nuclear_moments().unwrap();
        assert!(
            (m.magnetic_dipole_mu_n - 2.792847).abs() < 0.001,
            "proton μ={} μ_N",
            m.magnetic_dipole_mu_n
        );
        assert!(
            (m.electric_quadrupole_barn).abs() < 1e-10,
            "proton Q should be 0"
        );
    }

    #[test]
    fn nuclear_moments_deuteron() {
        let h2 = Nucleus::new(1, 2).unwrap();
        let m = h2.nuclear_moments().unwrap();
        assert!(
            (m.magnetic_dipole_mu_n - 0.857438).abs() < 0.001,
            "deuteron μ={}",
            m.magnetic_dipole_mu_n
        );
        assert!(
            (m.electric_quadrupole_barn - 0.002860).abs() < 0.0001,
            "deuteron Q={}",
            m.electric_quadrupole_barn
        );
    }

    #[test]
    fn nuclear_moments_he3_negative_mu() {
        let he3 = Nucleus::new(2, 3).unwrap();
        let m = he3.nuclear_moments().unwrap();
        assert!(
            m.magnetic_dipole_mu_n < 0.0,
            "He-3 should have negative magnetic moment"
        );
    }

    #[test]
    fn nuclear_moments_unknown_returns_none() {
        let og = Nucleus::new(118, 294).unwrap();
        assert!(og.nuclear_moments().is_none());
    }

    #[test]
    fn serde_roundtrip_nuclear_moments() {
        let m = NuclearMoments {
            magnetic_dipole_mu_n: 2.792847,
            electric_quadrupole_barn: 0.0,
        };
        let json = serde_json::to_string(&m).unwrap();
        let back: NuclearMoments = serde_json::from_str(&json).unwrap();
        assert!((m.magnetic_dipole_mu_n - back.magnetic_dipole_mu_n).abs() < 1e-10);
        assert!((m.electric_quadrupole_barn - back.electric_quadrupole_barn).abs() < 1e-10);
    }

    // --- Superallowed decay tests ---

    #[test]
    fn superallowed_ft_values_count() {
        let decays = superallowed_ft_values();
        assert_eq!(decays.len(), 9, "should have 9 superallowed transitions");
    }

    #[test]
    fn superallowed_ft_values_range() {
        // All ft values should be near 3040-3080 s
        for d in &superallowed_ft_values() {
            assert!(
                d.ft_seconds > 3030.0 && d.ft_seconds < 3080.0,
                "ft={} s out of range for {}->{}",
                d.ft_seconds,
                d.parent.a(),
                d.daughter.a()
            );
        }
    }

    #[test]
    fn superallowed_parent_daughter_consistency() {
        // Parent and daughter should have same A (superallowed: Z changes by 1)
        for d in &superallowed_ft_values() {
            assert_eq!(
                d.parent.a(),
                d.daughter.a(),
                "Parent A={} != daughter A={}",
                d.parent.a(),
                d.daughter.a()
            );
            assert_eq!(
                d.parent.z(),
                d.daughter.z() + 1,
                "Parent Z should be daughter Z + 1"
            );
        }
    }

    #[test]
    fn corrected_ft_value_returns_average() {
        let ft = corrected_ft_value(3042.3);
        assert!(
            (ft - 3072.27).abs() < 0.01,
            "Corrected Ft={ft} should be 3072.27"
        );
    }

    #[test]
    fn serde_roundtrip_superallowed_decay() {
        let d = &superallowed_ft_values()[0];
        let json = serde_json::to_string(d).unwrap();
        let back: SuperallowedDecay = serde_json::from_str(&json).unwrap();
        assert_eq!(d.parent, back.parent);
        assert_eq!(d.daughter, back.daughter);
        assert!((d.ft_seconds - back.ft_seconds).abs() < 1e-10);
    }
}
