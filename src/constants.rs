//! Physical constants from CODATA 2022.
//!
//! All values are taken from the 2022 CODATA recommended values of the
//! fundamental physical constants, published by NIST.
//!
//! Reference: <https://physics.nist.gov/cuu/Constants/>

/// Electron rest mass in MeV/c^2 (CODATA 2022).
pub const ELECTRON_MASS_MEV: f64 = 0.510_998_950_69;

/// Proton rest mass in MeV/c^2 (CODATA 2022).
pub const PROTON_MASS_MEV: f64 = 938.272_089_43;

/// Neutron rest mass in MeV/c^2 (CODATA 2022).
pub const NEUTRON_MASS_MEV: f64 = 939.565_421_94;

/// Atomic mass unit in MeV/c^2 (CODATA 2022).
pub const AMU_MEV: f64 = 931.494_103_72;

/// Fine-structure constant (CODATA 2022).
///
/// alpha = 1 / 137.035999177
pub const FINE_STRUCTURE: f64 = 1.0 / 137.035_999_177;

/// Reduced Planck constant in eV*s (CODATA 2022).
pub const HBAR_EV_S: f64 = 6.582_119_569e-16;

/// Rydberg constant in m^-1 (CODATA 2022).
pub const RYDBERG: f64 = 1.097_373_156_815_7e7;

/// Bohr radius in meters (CODATA 2022).
pub const BOHR_RADIUS: f64 = 5.291_772_105_44e-11;

/// Avogadro constant in mol^-1 (CODATA 2022, exact).
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// Elementary charge in coulombs (CODATA 2022, exact).
pub const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19;

/// Speed of light in vacuum in m/s (exact).
pub const C: f64 = 299_792_458.0;

/// Nuclear radius parameter r0 in femtometers.
///
/// R = r0 * A^(1/3), with r0 approximately 1.2 fm.
pub const R0_FM: f64 = 1.2;

/// Boltzmann constant in eV/K (CODATA 2022, exact).
pub const BOLTZMANN_EV: f64 = 8.617_333_262e-5;

/// Planck constant in eV*s (CODATA 2022, exact).
pub const H_EV_S: f64 = 4.135_667_696e-15;

/// Reduced Planck constant in MeV*s (CODATA 2022, exact).
///
/// Used for decay width ↔ lifetime conversion: Γ·τ = ħ.
pub const HBAR_MEV_S: f64 = 6.582_119_569e-22;

/// Bohr magneton in eV/T (CODATA 2022).
pub const BOHR_MAGNETON_EV_T: f64 = 5.788_381_798_2e-5;

/// Coulomb constant e²/(4πε₀) in MeV*fm.
///
/// e²/(4πε₀) = α ħc = 1.4399645 MeV*fm (CODATA 2022).
pub const COULOMB_MEV_FM: f64 = 1.439_964_546;

/// Classical electron radius in femtometers (CODATA 2022).
///
/// r_e = α ħc / (m_e c²) = 2.8179403205 fm.
///
/// Reference: <https://physics.nist.gov/cgi-bin/cuu/Value?re>
pub const CLASSICAL_ELECTRON_RADIUS_FM: f64 = 2.817_940_320_5;

/// Proton magnetic g-factor (CODATA 2022).
///
/// g_p = 2 μ_p / μ_N = 5.5856946893 (dimensionless).
///
/// Reference: <https://physics.nist.gov/cgi-bin/cuu/Value?gp>
pub const PROTON_G_FACTOR: f64 = 5.585_694_689_3;

/// Electron anomalous magnetic moment (CODATA 2022).
///
/// a_e = (|g_e| - 2) / 2 = 0.00115965218046.
///
/// Reference: <https://physics.nist.gov/cgi-bin/cuu/Value?ae>
pub const ELECTRON_ANOMALOUS_MOMENT: f64 = 0.001_159_652_180_46;

/// Standard acceleration of gravity in m/s² (CODATA 2022, exact).
///
/// Reference: <https://physics.nist.gov/cgi-bin/cuu/Value?gn>
pub const STANDARD_GRAVITY: f64 = 9.806_65;

/// Geocentric gravitational constant GM_earth in m³/s² (IERS 2010).
///
/// Reference: IERS Conventions (2010), Table 1.1.
pub const GM_EARTH: f64 = 3.986_004_418e14;

/// Earth angular rotation rate in rad/s (IERS Conventions 2010, Table 1.2).
///
/// Reference: IERS Conventions (2010), Table 1.2.
pub const EARTH_ROTATION_RAD_S: f64 = 7.292_115_0e-5;

/// Reduced Planck constant times c in MeV*fm (CODATA 2022, exact).
pub const HBAR_C_MEV_FM: f64 = 197.326_980_4;

/// Planck constant times c in eV*nm (CODATA 2022, exact).
pub const HC_EV_NM: f64 = 1_239.841_984;

/// Rydberg energy R∞hc in eV (CODATA 2022).
pub const RYDBERG_EV: f64 = 13.605_693_122_990;

/// Electron mass in atomic mass units (CODATA 2022).
pub const ELECTRON_MASS_U: f64 = 5.485_799_090_441e-4;

/// Proton mass in atomic mass units (CODATA 2022).
pub const PROTON_MASS_U: f64 = 1.007_276_466_578_9;

/// Neutron mass in atomic mass units (CODATA 2022).
pub const NEUTRON_MASS_U: f64 = 1.008_664_916_06;

/// Hydrogen-1 atomic mass in atomic mass units (AME2020).
pub const HYDROGEN_ATOM_MASS_U: f64 = 1.007_825_031_898;

/// Nuclear magneton in eV/T (CODATA 2022).
pub const NUCLEAR_MAGNETON_EV_T: f64 = 3.152_451_254_17e-8;

/// Deuteron magnetic g-factor g_d = μ_d/μ_N for spin 1 (CODATA 2022).
pub const DEUTERON_G_FACTOR: f64 = 0.857_438_233_5;

/// Atomic unit of time ħ/E_h in seconds (CODATA 2022).
pub const ATOMIC_UNIT_TIME_S: f64 = 2.418_884_326_586_4e-17;

/// Coulomb constant 1/(4πε₀) in N*m²/C² (CODATA 2022).
pub const COULOMB_K_SI: f64 = 8.987_551_786_2e9;
