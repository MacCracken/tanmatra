//! Physical constants from CODATA 2022.
//!
//! All values are taken from the 2022 CODATA recommended values of the
//! fundamental physical constants, published by NIST.
//!
//! Reference: <https://physics.nist.gov/cuu/Constants/>

/// Electron rest mass in MeV/c^2 (CODATA 2022).
pub const ELECTRON_MASS_MEV: f64 = 0.510_998_950;

/// Proton rest mass in MeV/c^2 (CODATA 2022).
pub const PROTON_MASS_MEV: f64 = 938.272_088_16;

/// Neutron rest mass in MeV/c^2 (CODATA 2022).
pub const NEUTRON_MASS_MEV: f64 = 939.565_420_52;

/// Atomic mass unit in MeV/c^2 (CODATA 2022).
pub const AMU_MEV: f64 = 931.494_102_42;

/// Fine-structure constant (CODATA 2022).
///
/// alpha = 1 / 137.035999084
pub const FINE_STRUCTURE: f64 = 1.0 / 137.035_999_084;

/// Reduced Planck constant in eV*s (CODATA 2022).
pub const HBAR_EV_S: f64 = 6.582_119_569e-16;

/// Rydberg constant in m^-1 (CODATA 2022).
pub const RYDBERG: f64 = 1.097_373_156_816_0e7;

/// Bohr radius in meters (CODATA 2022).
pub const BOHR_RADIUS: f64 = 5.291_772_109_03e-11;

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

/// Reduced Planck constant in MeV*s (CODATA 2022).
///
/// Used for decay width ↔ lifetime conversion: Γ·τ = ħ.
pub const HBAR_MEV_S: f64 = 6.582_119_514e-22;

/// Bohr magneton in eV/T (CODATA 2022).
pub const BOHR_MAGNETON_EV_T: f64 = 5.788_381_806_0e-5;

/// Coulomb constant k_e in MeV*fm / e^2.
///
/// k_e = e^2 / (4 * pi * epsilon_0) expressed in nuclear units.
/// Numerically: k_e * e^2 = 1.4399764 MeV*fm.
pub const COULOMB_MEV_FM: f64 = 1.439_964_5;
