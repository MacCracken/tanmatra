//! Nuclear structure and binding energy calculations.
//!
//! Implements the semi-empirical mass formula (Bethe-Weizsacker formula) for
//! nuclear binding energies, along with nuclear radius calculations and
//! magic number identification.

use crate::constants::{AMU_MEV, NEUTRON_MASS_MEV, PROTON_MASS_MEV, R0_FM};
use crate::error::TanmatraError;
use serde::{Deserialize, Serialize};

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

/// Returns `true` if the given number is a nuclear magic number.
///
/// Magic numbers correspond to complete nuclear shells:
/// 2, 8, 20, 28, 50, 82, 126.
#[must_use]
pub fn is_magic_number(n: u32) -> bool {
    matches!(n, 2 | 8 | 20 | 28 | 50 | 82 | 126)
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
}
