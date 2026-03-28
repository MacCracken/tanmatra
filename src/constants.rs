//! CODATA 2022 fundamental physical constants.
//!
//! All values sourced from NIST CODATA 2022 recommended values.

/// Electron rest mass in MeV/c^2.
pub const ELECTRON_MASS_MEV: f64 = 0.510_998_95;

/// Proton rest mass in MeV/c^2.
pub const PROTON_MASS_MEV: f64 = 938.272_088_16;

/// Neutron rest mass in MeV/c^2.
pub const NEUTRON_MASS_MEV: f64 = 939.565_420_52;

/// Atomic mass unit in MeV/c^2.
pub const ATOMIC_MASS_UNIT_MEV: f64 = 931.494_102_42;

/// Fine-structure constant (dimensionless).
pub const FINE_STRUCTURE: f64 = 1.0 / 137.035_999_084;

/// Reduced Planck constant in eV*s.
pub const HBAR_EV_S: f64 = 6.582_119_569e-16;

/// Rydberg constant in 1/m.
pub const RYDBERG_CONST: f64 = 1.097_373_156_816_0e7;

/// Bohr radius in meters.
pub const BOHR_RADIUS_M: f64 = 5.291_772_109_03e-11;

/// Avogadro constant in 1/mol.
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// Elementary charge in coulombs.
pub const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19;

/// Speed of light in vacuum in m/s.
pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;

/// Boltzmann constant in J/K.
pub const BOLTZMANN: f64 = 1.380_649e-23;

/// Nuclear radius constant in femtometers.
/// R = R0 * A^(1/3), where R0 is approximately 1.2 fm.
pub const R0_FM: f64 = 1.2;

/// Coulomb constant `k_e` in `MeV`*fm / e^2.
/// Numerically: `1.439_964_5` `MeV`*fm.
pub const COULOMB_CONST_MEV_FM: f64 = 1.439_964_5;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fine_structure_approx() {
        let expected = 1.0 / 137.036;
        let diff = libm::fabs(FINE_STRUCTURE - expected);
        assert!(diff < 1e-6);
    }

    #[test]
    fn proton_heavier_than_electron() {
        const { assert!(PROTON_MASS_MEV > ELECTRON_MASS_MEV) };
    }

    #[test]
    fn neutron_heavier_than_proton() {
        const { assert!(NEUTRON_MASS_MEV > PROTON_MASS_MEV) };
    }

    #[test]
    fn speed_of_light_exact() {
        assert!((SPEED_OF_LIGHT - 299_792_458.0).abs() < f64::EPSILON);
    }
}
