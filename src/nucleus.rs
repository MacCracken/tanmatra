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
    (1, 1, 7288.971064),    // H-1
    (1, 2, 13135.722895),   // H-2
    (2, 3, 14931.21888),    // He-3
    (2, 4, 2424.91587),     // He-4
    (3, 6, 14086.88044),    // Li-6
    (3, 7, 14907.10463),    // Li-7
    (6, 12, 0.0),           // C-12 (definition)
    (6, 13, 3125.00933),    // C-13
    (7, 14, 2863.41683),    // N-14
    (8, 16, -4737.00217),   // O-16
    (9, 19, -1487.44512),   // F-19
    (11, 23, -9529.85352),  // Na-23
    (14, 28, -21492.79711), // Si-28
    (15, 31, -24440.54442), // P-31
    (16, 32, -26015.53714), // S-32
    (20, 40, -34846.402),   // Ca-40
    (26, 56, -60607.163),   // Fe-56
    (28, 58, -60228.871),   // Ni-58
    (28, 62, -66746.440),   // Ni-62
    (29, 63, -65579.868),   // Cu-63
    (30, 64, -66004.018),   // Zn-64
    (38, 88, -87921.62876), // Sr-88
    (40, 90, -88772.547),   // Zr-90
    (42, 98, -88115.980),   // Mo-98
    (50, 120, -91097.741),  // Sn-120
    (53, 127, -88983.217),  // I-127
    (55, 133, -88070.943),  // Cs-133
    (56, 138, -88261.806),  // Ba-138
    (82, 208, -21748.519),  // Pb-208
    (83, 209, -18258.589),  // Bi-209
    (92, 235, 40918.782),   // U-235
    (92, 238, 47307.732),   // U-238
];

/// Conversion factor used by AME2020: 1 u = 931494.10242 keV (CODATA 2018,
/// the value AME2020 was evaluated with).
const AME2020_U_KEV: f64 = 931_494.102_42;

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
    (28, 58, 3.7757),  // Ni-58
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
/// Sources: magnetic dipole moments — N.J. Stone, "Table of recommended
/// nuclear magnetic dipole moments", INDC(NDS)-0794 (IAEA, 2019); electric
/// quadrupole moments — N.J. Stone, "Table of recommended nuclear electric
/// quadrupole moments", INDC(NDS)-0833 (IAEA, 2021). Nuclei with ground-state
/// spin 0 or 1/2 have no spectroscopic quadrupole moment (Q = 0).
const NUCLEAR_MOMENTS: &[(u32, u32, f64, f64)] = &[
    (1, 1, 2.792847351, 0.0),       // H-1   (1/2+)
    (1, 2, 0.857438231, 0.0028578), // H-2   (1+)
    (2, 3, -2.12762531, 0.0),       // He-3  (1/2+)
    (3, 6, 0.822043, -0.000806),    // Li-6  (1+)
    (3, 7, 3.256407, -0.0400),      // Li-7  (3/2-)
    (6, 13, 0.702369, 0.0),         // C-13  (1/2-)
    (7, 14, 0.403573, 0.02044),     // N-14  (1+)
    (8, 17, -1.893543, -0.0256),    // O-17  (5/2+)
    (9, 19, 2.628321, 0.0),         // F-19  (1/2+: no spectroscopic Q)
    (11, 23, 2.21750, 0.104),       // Na-23 (3/2+)
    (13, 27, 3.64070, 0.1466),      // Al-27 (5/2+)
    (15, 31, 1.130925, 0.0),        // P-31  (1/2+)
    (55, 133, 2.5778, -0.00343),    // Cs-133 (7/2+)
    (82, 207, 0.5906, 0.0),         // Pb-207 (1/2-)
    (83, 209, 4.092, -0.516),       // Bi-209 (9/2-)
    (92, 235, -0.38, 4.936),        // U-235 (7/2-)
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

impl SuperallowedDecay {
    /// Returns the corrected Ft = ft(1 + δ'_R)(1 + δ_NS − δ_C) in seconds for
    /// this transition, from Hardy & Towner 2020 Table XVI.
    ///
    /// Returns `None` if the transition is not in the table.
    #[must_use]
    pub fn corrected_ft_seconds(&self) -> Option<f64> {
        SUPERALLOWED_FT_VALUES
            .iter()
            .find(|&&(pz, pa, dz, da, _, _)| {
                pz == self.parent.z()
                    && pa == self.parent.a()
                    && dz == self.daughter.z()
                    && da == self.daughter.a()
            })
            .map(|&(_, _, _, _, _, big_ft)| big_ft)
    }
}

/// Superallowed 0+ → 0+ beta-decay ft values.
///
/// Source: Hardy, J.C. & Towner, I.S., Physical Review C 102, 045501 (2020),
/// Table XVI: measured ft and corrected Ft = ft(1 + δ'_R)(1 + δ_NS − δ_C).
/// The Al-26 and K-38 parents are the 0+ isomers; `Nucleus` identifies only (Z, A).
const SUPERALLOWED_FT_VALUES: &[(u32, u32, u32, u32, f64, f64)] = &[
    // (parent_z, parent_a, daughter_z, daughter_a, ft_seconds, corrected_Ft_seconds)
    (8, 14, 7, 14, 3042.2, 3070.2), // O-14 -> N-14 (2.313 MeV 0+ state)
    (13, 26, 12, 26, 3037.61, 3072.4), // Al-26m (228 keV isomer) -> Mg-26
    (17, 34, 16, 34, 3049.43, 3071.6), // Cl-34 -> S-34
    (19, 38, 18, 38, 3051.45, 3072.9), // K-38m (130 keV isomer) -> Ar-38
    (21, 42, 20, 42, 3047.7, 3071.7), // Sc-42 -> Ca-42
    (23, 46, 22, 46, 3050.33, 3074.3), // V-46 -> Ti-46
    (25, 50, 24, 50, 3048.4, 3071.1), // Mn-50 -> Cr-50
    (27, 54, 26, 54, 3050.8, 3070.4), // Co-54 -> Fe-54
    (31, 62, 30, 62, 3074.1, 3072.4), // Ga-62 -> Zn-62
];

/// World-average corrected Ft from Hardy & Towner 2020, eq. (22):
/// Ft = 3072.24 ± 0.57(stat) ± 0.36(δ_NS) → ± 1.85 s total.
///
/// Source: Hardy, J.C. & Towner, I.S., Physical Review C 102, 045501 (2020).
const AVERAGE_CORRECTED_FT: f64 = 3072.24;

/// Bethe-Weizsacker semi-empirical mass formula coefficients (in MeV).
///
/// Least-squares fit of B = a_v A − a_s A^(2/3) − a_c Z(Z−1)/A^(1/3)
/// − a_a (A−2Z)²/A + a_p δ/√A to all 2484 experimental (non-systematic)
/// AME2020 binding energies with A ≥ 16. RMS deviation 3.31 MeV.
const A_V: f64 = 15.413_751; // Volume term
const A_S: f64 = 16.860_179; // Surface term
const A_C: f64 = 0.695_241; // Coulomb term
const A_A: f64 = 22.497_475; // Asymmetry term
/// Pairing term coefficient (MeV).
const A_P: f64 = 12.027_823;

/// Myers–Swiatecki shell-correction amplitude C (MeV), fitted to the AME2020
/// residuals of the liquid-drop fit above with the published c = 0.325.
/// RMS deviation with the correction: 2.76 MeV.
const MS_SHELL_C: f64 = 2.886_704;
/// Myers–Swiatecki shell-correction constant c (Myers & Swiatecki,
/// Nucl. Phys. 81, 1 (1966)).
const MS_SHELL_SMALL_C: f64 = 0.325;

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

    /// Calculates binding energy with the Myers–Swiatecki shell correction.
    ///
    /// B_corrected = B_LDM − S(N, Z), with the shell term of Myers & Swiatecki,
    /// Nucl. Phys. 81, 1 (1966):
    ///
    /// S = C [ (F(N) + F(Z)) / (A/2)^(2/3) − c A^(1/3) ]
    ///
    /// F(N) = (3/5) [(M_i^(5/3) − M_(i−1)^(5/3)) / (M_i − M_(i−1))] (N − M_(i−1))
    ///        − (3/5) (N^(5/3) − M_(i−1)^(5/3)),  for M_(i−1) < N ≤ M_i,
    ///
    /// with magic numbers M = 0, 2, 8, 20, 28, 50, 82, 126, 184. S is negative
    /// (extra binding) at closed shells and positive between them. C is fitted
    /// to AME2020 (2.8867 MeV); c = 0.325 as published. Against AME2020 (A ≥ 16)
    /// the RMS deviation drops from 3.31 MeV to 2.76 MeV.
    #[must_use]
    #[inline]
    pub fn binding_energy_shell_corrected(&self) -> f64 {
        if self.a == 1 {
            return 0.0;
        }
        self.binding_energy() - myers_swiatecki_shell_term(self.z, self.n())
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

    /// Returns the atomic mass in atomic mass units (u) from the
    /// semi-empirical nuclear mass.
    ///
    /// M_atom = M_nucleus + Z m_e − B_e(Z), where the total electron binding
    /// energy is B_e(Z) = 14.4381 Z^2.39 + 1.55468e-6 Z^5.35 eV
    /// (Lunney, Pearson & Thibault, Rev. Mod. Phys. 75, 1021 (2003), eq. A4).
    #[must_use]
    #[inline]
    pub fn atomic_mass_amu(&self) -> f64 {
        let z = self.z as f64;
        let electron_binding_mev =
            (14.4381 * libm::pow(z, 2.39) + 1.554_68e-6 * libm::pow(z, 5.35)) * 1e-6;
        (self.nuclear_mass() + z * crate::constants::ELECTRON_MASS_MEV - electron_binding_mev)
            / AMU_MEV
    }

    /// Returns the nuclear mass in atomic mass units (u).
    #[must_use]
    #[inline]
    pub fn nuclear_mass_amu(&self) -> f64 {
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

/// Magic numbers used by the Myers–Swiatecki shell function.
const MS_MAGIC: [u32; 9] = [0, 2, 8, 20, 28, 50, 82, 126, 184];

/// x^(5/3) computed as x·cbrt(x)².
fn pow_five_thirds(x: f64) -> f64 {
    let c = libm::cbrt(x);
    x * c * c
}

/// Myers–Swiatecki F(N) for one nucleon species (dimensionless).
fn myers_swiatecki_f(count: u32) -> f64 {
    for w in MS_MAGIC.windows(2) {
        let (lo, hi) = (w[0], w[1]);
        if count <= hi {
            let lo_f = lo as f64;
            let hi_f = hi as f64;
            let n = count as f64;
            let lo_53 = pow_five_thirds(lo_f);
            let q = 0.6 * (pow_five_thirds(hi_f) - lo_53) / (hi_f - lo_f);
            return q * (n - lo_f) - 0.6 * (pow_five_thirds(n) - lo_53);
        }
    }
    0.0
}

/// Myers–Swiatecki shell term S(N, Z) in MeV (positive = less bound).
fn myers_swiatecki_shell_term(z: u32, n: u32) -> f64 {
    let a = (z + n) as f64;
    let c = libm::cbrt(a / 2.0);
    let half_a_two_thirds = c * c;
    MS_SHELL_C
        * ((myers_swiatecki_f(n) + myers_swiatecki_f(z)) / half_a_two_thirds
            - MS_SHELL_SMALL_C * libm::cbrt(a))
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

/// Empirical single-particle filling order for protons.
///
/// Same levels and closures as `SHELL_MODEL_LEVELS`, with the order inside
/// each major shell chosen to maximise agreement with the firm ground-state
/// J^π of odd-Z, even-N nuclei in NUBASE2020 (243 of 427).
const PROTON_LEVEL_ORDER: [(u32, u32, u32); 32] = [
    (1, 0, 1),
    (1, 1, 3),
    (1, 1, 1),
    (1, 2, 5),
    (2, 0, 1),
    (1, 2, 3),
    (1, 3, 7),
    (2, 1, 3),
    (1, 3, 5),
    (2, 1, 1),
    (1, 4, 9),
    (1, 4, 7),
    (2, 2, 5),
    (1, 5, 11),
    (2, 2, 3),
    (3, 0, 1),
    (1, 5, 9),
    (2, 3, 5),
    (2, 3, 7),
    (3, 1, 3),
    (3, 1, 1),
    (1, 6, 13),
    (2, 4, 9),
    (1, 6, 11),
    (3, 2, 5),
    (4, 0, 1),
    (2, 4, 7),
    (3, 2, 3),
    (1, 7, 15),
    (2, 5, 11),
    (2, 5, 9),
    (3, 3, 7),
];

/// Empirical single-particle filling order for neutrons.
///
/// Chosen as for `PROTON_LEVEL_ORDER` against the odd-N, even-Z nuclei of
/// NUBASE2020 (172 of 421). Reproduces e.g. Pb-207 = 1/2⁻ (3p1/2 hole).
const NEUTRON_LEVEL_ORDER: [(u32, u32, u32); 32] = [
    (1, 0, 1),
    (1, 1, 3),
    (1, 1, 1),
    (1, 2, 5),
    (2, 0, 1),
    (1, 2, 3),
    (1, 3, 7),
    (2, 1, 3),
    (1, 3, 5),
    (2, 1, 1),
    (1, 4, 9),
    (2, 2, 5),
    (3, 0, 1),
    (1, 4, 7),
    (1, 5, 11),
    (2, 2, 3),
    (2, 3, 7),
    (3, 1, 3),
    (1, 6, 13),
    (1, 5, 9),
    (2, 3, 5),
    (3, 1, 1),
    (2, 4, 9),
    (1, 6, 11),
    (3, 2, 5),
    (4, 0, 1),
    (2, 4, 7),
    (3, 2, 3),
    (1, 7, 15),
    (2, 5, 11),
    (2, 5, 9),
    (3, 3, 7),
];

/// Returns the level holding the last nucleon for `count` nucleons filled in
/// `order`, and how many nucleons occupy it.
fn last_filled_level(order: &[(u32, u32, u32)], count: u32) -> Option<(ShellLevel, u32)> {
    let mut remaining = count;
    for &(n_shell, l, two_j) in order {
        let level = ShellLevel { n_shell, l, two_j };
        let deg = level.degeneracy();
        if remaining <= deg {
            return Some((level, remaining));
        }
        remaining -= deg;
    }
    None
}

/// Returns the ground-state spin and parity (J^pi) of a nucleus in the
/// extreme single-particle shell model.
///
/// Returns `(two_j, parity)` where `parity` is +1 or −1 and `two_j` is twice
/// the nuclear spin (even for even A, odd for odd A).
///
/// - Even-even: 0⁺.
/// - Odd A: j and parity (−1)^l of the level holding the unpaired nucleon,
///   using empirical proton and neutron level orders (the Mayer–Jensen levels
///   of [`shell_model_levels`], reordered within each major shell to best match
///   NUBASE2020 ground states).
/// - Odd-odd: parity (−1)^(l_p + l_n); spin from the Brennan–Bernstein
///   coupling rules (Phys. Rev. 120, 927 (1960), extending Nordheim,
///   Rev. Mod. Phys. 23, 322 (1951)). A particle–hole pair (one level less
///   than half full, the other more) gives J = j_p + j_n − 1. Otherwise, with
///   Nordheim number N = (j_p − l_p) + (j_n − l_n): N = 0 gives
///   J = |j_p − j_n|, N = ±1 gives J = j_p + j_n.
///
/// Agreement with firm NUBASE2020 ground states: odd A 415/848, odd-odd
/// 107/352 (parity 267/352). The extreme single-particle model does not
/// describe deformed nuclei.
///
/// Nucleon counts above 184 use the last listed level.
#[must_use]
pub fn ground_state_spin_parity(nucleus: &Nucleus) -> (u32, i32) {
    let z = nucleus.z();
    let n = nucleus.n();
    let z_odd = !z.is_multiple_of(2);
    let n_odd = !n.is_multiple_of(2);
    let parity_of = |l: u32| if l.is_multiple_of(2) { 1 } else { -1 };
    let last = |order: &[(u32, u32, u32)], count: u32| {
        last_filled_level(order, count).unwrap_or_else(|| {
            let (n_shell, l, two_j) = order[order.len() - 1];
            (ShellLevel { n_shell, l, two_j }, 1)
        })
    };

    match (z_odd, n_odd) {
        (false, false) => (0, 1),
        (true, false) => {
            let (p, _) = last(&PROTON_LEVEL_ORDER, z);
            (p.two_j, parity_of(p.l))
        }
        (false, true) => {
            let (nl, _) = last(&NEUTRON_LEVEL_ORDER, n);
            (nl.two_j, parity_of(nl.l))
        }
        (true, true) => {
            let (p, p_occ) = last(&PROTON_LEVEL_ORDER, z);
            let (nl, n_occ) = last(&NEUTRON_LEVEL_ORDER, n);
            let parity = parity_of(p.l + nl.l);
            let p_hole = p_occ > p.degeneracy() / 2;
            let n_hole = n_occ > nl.degeneracy() / 2;
            // 2(j − l) = +1 or −1.
            let sign_p = p.two_j as i64 - 2 * p.l as i64;
            let sign_n = nl.two_j as i64 - 2 * nl.l as i64;
            let two_j = if p_hole != n_hole {
                p.two_j + nl.two_j - 2
            } else if sign_p + sign_n == 0 {
                (p.two_j as i64 - nl.two_j as i64).unsigned_abs() as u32
            } else {
                p.two_j + nl.two_j
            };
            (two_j, parity)
        }
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
    for &(pz, pa, dz, da, ft, _) in SUPERALLOWED_FT_VALUES {
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

/// Returns the world-average corrected Ft value, 3072.24 s (Hardy & Towner 2020).
///
/// The corrections δ'_R, δ_NS and δ_C are transition-specific and cannot be
/// derived from a bare ft value, so the argument is not used. For the
/// corrected Ft of a specific transition use
/// [`SuperallowedDecay::corrected_ft_seconds`].
///
/// Source: Hardy & Towner, Phys. Rev. C 102, 045501 (2020), eq. (22).
#[must_use]
#[deprecated(
    since = "1.3.0",
    note = "the argument is ignored; use SuperallowedDecay::corrected_ft_seconds or superallowed_average_ft"
)]
pub fn corrected_ft_value(_ft: f64) -> f64 {
    AVERAGE_CORRECTED_FT
}

/// Returns the world-average corrected Ft value in seconds and its total
/// uncertainty: (3072.24, 1.85) (Hardy & Towner 2020, eq. 22).
#[must_use]
pub const fn superallowed_average_ft() -> (f64, f64) {
    (AVERAGE_CORRECTED_FT, 1.85)
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

    // --- Shell correction tests ---

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
        assert!(
            (me - (-60607.163)).abs() < 1e-6,
            "Fe-56 mass excess={me} keV"
        );
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
    #[allow(deprecated)]
    fn corrected_ft_value_returns_average() {
        let ft = corrected_ft_value(3042.3);
        assert!((ft - 3072.24).abs() < 1e-9, "Corrected Ft={ft}");
        assert_eq!(superallowed_average_ft(), (3072.24, 1.85));
    }

    #[test]
    fn superallowed_corrected_ft_per_transition() {
        let decays = superallowed_ft_values();
        let o14 = decays.iter().find(|d| d.parent.a() == 14).unwrap();
        assert!((o14.ft_seconds - 3042.2).abs() < 1e-9);
        assert!((o14.corrected_ft_seconds().unwrap() - 3070.2).abs() < 1e-9);
        let ga62 = decays.iter().find(|d| d.parent.a() == 62).unwrap();
        assert!((ga62.corrected_ft_seconds().unwrap() - 3072.4).abs() < 1e-9);
    }

    // --- Reference-value tests against AME2020 / NUBASE2020 ---

    #[test]
    fn semf_close_to_ame2020() {
        // AME2020 binding energies (MeV): B/A x A.
        for (z, a, b_exp) in [
            (26, 56, 492.2600),
            (82, 208, 1636.4301),
            (50, 120, 1020.5448),
        ] {
            let b = Nucleus::new(z, a).unwrap().binding_energy();
            assert!((b - b_exp).abs() < 12.0, "Z={z} A={a}: B={b} vs {b_exp}");
        }
    }

    #[test]
    fn shell_correction_improves_doubly_magic() {
        // AME2020: Ni-56, Sn-132, Pb-208.
        for (z, a, b_exp) in [
            (28, 56, 483.9956),
            (50, 132, 1102.8432),
            (82, 208, 1636.4301),
        ] {
            let nuc = Nucleus::new(z, a).unwrap();
            let err_ldm = (nuc.binding_energy() - b_exp).abs();
            let err_shell = (nuc.binding_energy_shell_corrected() - b_exp).abs();
            assert!(err_shell < err_ldm, "Z={z} A={a}: {err_shell} vs {err_ldm}");
        }
    }

    #[test]
    fn shell_term_sign_at_magic_and_midshell() {
        // Doubly magic Pb-208: extra binding. Mid-shell Er-166: less binding.
        assert!(myers_swiatecki_shell_term(82, 126) < 0.0);
        assert!(myers_swiatecki_shell_term(68, 98) > 0.0);
    }

    #[test]
    fn atomic_mass_includes_electrons() {
        let fe = Nucleus::iron_56();
        let diff = fe.atomic_mass_amu() - fe.nuclear_mass_amu();
        // 26 electrons (0.0142631 u) minus the 34.8 keV total electron
        // binding energy (0.0000374 u).
        assert!((diff - 0.014_225_7).abs() < 1e-7, "diff={diff}");
    }

    #[test]
    fn odd_odd_spin_is_integer() {
        for z in 1..=99u32 {
            for a in (2 * z)..=(2 * z + 60) {
                if z % 2 == 1 && (a - z) % 2 == 1 {
                    let (two_j, _) = ground_state_spin_parity(&Nucleus::new(z, a).unwrap());
                    assert_eq!(two_j % 2, 0, "Z={z} A={a} gave half-integer spin");
                }
            }
        }
    }

    #[test]
    fn spin_parity_reference_nuclei() {
        // NUBASE2020 ground states.
        for (z, a, two_j, parity) in [
            (7, 14, 2, 1),    // N-14 1+
            (5, 10, 6, 1),    // B-10 3+
            (19, 40, 8, -1),  // K-40 4-
            (17, 38, 4, -1),  // Cl-38 2-
            (82, 207, 1, -1), // Pb-207 1/2-
            (83, 209, 9, -1), // Bi-209 9/2-
            (8, 17, 5, 1),    // O-17 5/2+
            (20, 41, 7, -1),  // Ca-41 7/2-
        ] {
            let got = ground_state_spin_parity(&Nucleus::new(z, a).unwrap());
            assert_eq!(got, (two_j, parity), "Z={z} A={a}");
        }
    }

    #[test]
    fn spin_zero_or_half_has_no_quadrupole() {
        for &(z, a, _, q) in NUCLEAR_MOMENTS {
            let (two_j, _) = match (z, a) {
                (9, 19) | (1, 1) | (2, 3) | (6, 13) | (15, 31) | (82, 207) => (1, 0),
                _ => (2, 0),
            };
            if two_j <= 1 {
                assert!(q.abs() < 1e-15, "Z={z} A={a} spin <= 1/2 but Q={q}");
            }
        }
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
