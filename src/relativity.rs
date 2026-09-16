//! Relativistic kinematics.
//!
//! Lorentz transformations, four-vectors, invariant mass, and
//! relativistic energy-momentum relations.

use crate::constants::C;
use serde::{Deserialize, Serialize};

/// A relativistic four-momentum (E/c, px, py, pz) in MeV/c units.
///
/// The invariant mass is: m²c⁴ = E² - (pc)² = E² - px² - py² - pz²
/// (when all components are in MeV).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct FourMomentum {
    /// Energy component E in MeV.
    pub energy: f64,
    /// x-component of momentum in MeV/c.
    pub px: f64,
    /// y-component of momentum in MeV/c.
    pub py: f64,
    /// z-component of momentum in MeV/c.
    pub pz: f64,
}

impl FourMomentum {
    /// Creates a four-momentum from energy and momentum components (all MeV).
    #[must_use]
    pub const fn new(energy: f64, px: f64, py: f64, pz: f64) -> Self {
        Self { energy, px, py, pz }
    }

    /// Creates a four-momentum for a particle at rest with given mass.
    #[must_use]
    pub const fn at_rest(mass_mev: f64) -> Self {
        Self {
            energy: mass_mev,
            px: 0.0,
            py: 0.0,
            pz: 0.0,
        }
    }

    /// Creates a four-momentum from mass and 3-momentum magnitude along z.
    #[must_use]
    pub fn from_mass_and_momentum(mass_mev: f64, p_mev: f64) -> Self {
        let energy = libm::sqrt(mass_mev * mass_mev + p_mev * p_mev);
        Self {
            energy,
            px: 0.0,
            py: 0.0,
            pz: p_mev,
        }
    }

    /// Returns the invariant mass in MeV/c².
    ///
    /// m = √(E² - p²) where all in MeV.
    #[must_use]
    #[inline]
    pub fn invariant_mass(&self) -> f64 {
        // (E − p)(E + p) avoids the cancellation in E² − p² when E ≫ m.
        let p = self.momentum_magnitude();
        let m_sq = (self.energy - p) * (self.energy + p);
        if m_sq < 0.0 {
            return 0.0; // spacelike (numerical noise for massless particles)
        }
        libm::sqrt(m_sq)
    }

    /// Returns the magnitude of the 3-momentum in MeV/c.
    #[must_use]
    #[inline]
    pub fn momentum_magnitude(&self) -> f64 {
        libm::hypot(libm::hypot(self.px, self.py), self.pz)
    }

    /// Returns the velocity β = v/c (dimensionless).
    #[must_use]
    #[inline]
    pub fn beta(&self) -> f64 {
        if self.energy <= 0.0 {
            return 0.0;
        }
        let p = self.momentum_magnitude();
        p / self.energy
    }

    /// Returns the Lorentz factor γ = E / (mc²).
    #[must_use]
    #[inline]
    pub fn gamma(&self) -> f64 {
        let mass = self.invariant_mass();
        if mass <= 0.0 {
            return f64::INFINITY; // massless particle
        }
        self.energy / mass
    }

    /// Returns the kinetic energy T = E - mc² in MeV.
    #[must_use]
    #[inline]
    pub fn kinetic_energy(&self) -> f64 {
        // T = E − m = p² / (E + m), without cancellation at low momentum.
        let m = self.invariant_mass();
        let p = self.momentum_magnitude();
        if self.energy + m <= 0.0 {
            return 0.0;
        }
        p * p / (self.energy + m)
    }
}

impl core::ops::Add for FourMomentum {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            energy: self.energy + rhs.energy,
            px: self.px + rhs.px,
            py: self.py + rhs.py,
            pz: self.pz + rhs.pz,
        }
    }
}

/// Calculates the Lorentz factor γ for a given velocity β = v/c.
///
/// γ = 1 / √(1 - β²)
///
/// Returns [`f64::INFINITY`] if |β| >= 1.
#[must_use]
#[inline]
pub fn lorentz_gamma(beta: f64) -> f64 {
    let b = libm::fabs(beta);
    if b >= 1.0 {
        return f64::INFINITY;
    }
    // (1 − β)(1 + β) keeps full precision as β → 1.
    1.0 / libm::sqrt((1.0 - b) * (1.0 + b))
}

/// Calculates the velocity β = v/c from the Lorentz factor γ.
///
/// β = √(1 - 1/γ²)
#[must_use]
#[inline]
pub fn gamma_to_beta(gamma: f64) -> f64 {
    if gamma <= 1.0 {
        return 0.0;
    }
    if gamma.is_infinite() {
        return 1.0;
    }
    libm::sqrt((gamma - 1.0) * (gamma + 1.0)) / gamma
}

/// Relativistic velocity addition for collinear velocities:
/// β_total = (β1 + β2) / (1 + β1*β2).
///
/// Both velocities must be in units of c with |β| ≤ 1. Returns NaN when the
/// result is undefined (|β| > 1, or light speeds in opposite directions,
/// where 1 + β1β2 = 0).
#[must_use]
#[inline]
pub fn velocity_addition(beta1: f64, beta2: f64) -> f64 {
    let denom = 1.0 + beta1 * beta2;
    if libm::fabs(beta1) > 1.0 || libm::fabs(beta2) > 1.0 || denom == 0.0 {
        return f64::NAN;
    }
    (beta1 + beta2) / denom
}

/// Calculates the relativistic total energy from mass and momentum.
///
/// E = √((mc²)² + (pc)²)
///
/// Parameters in MeV and MeV/c respectively. Returns MeV.
#[must_use]
#[inline]
pub fn relativistic_energy(mass_mev: f64, momentum_mev: f64) -> f64 {
    libm::sqrt(mass_mev * mass_mev + momentum_mev * momentum_mev)
}

/// Calculates the relativistic momentum from mass and kinetic energy.
///
/// p = √(T² + 2Tmc²) / c
///
/// Parameters in MeV. Returns MeV/c.
#[must_use]
#[inline]
pub fn relativistic_momentum(mass_mev: f64, kinetic_energy_mev: f64) -> f64 {
    libm::sqrt(kinetic_energy_mev * kinetic_energy_mev + 2.0 * kinetic_energy_mev * mass_mev)
}

/// Calculates the de Broglie wavelength in femtometers.
///
/// λ = h / p = hc / (pc)
///
/// Parameter: momentum in MeV/c. Returns wavelength in fm.
#[must_use]
#[inline]
pub fn de_broglie_wavelength_fm(momentum_mev: f64) -> f64 {
    if momentum_mev <= 0.0 {
        return f64::INFINITY;
    }
    // λ = 2π ħc / (pc)
    2.0 * core::f64::consts::PI * crate::constants::HBAR_C_MEV_FM / momentum_mev
}

/// Converts a velocity in m/s to β = v/c.
#[must_use]
#[inline]
pub fn velocity_to_beta(velocity_m_per_s: f64) -> f64 {
    velocity_m_per_s / C
}

/// Converts β = v/c to velocity in m/s.
#[must_use]
#[inline]
pub fn beta_to_velocity(beta: f64) -> f64 {
    beta * C
}

/// Calculates the invariant mass of a two-body system from their four-momenta.
///
/// M² = (p1 + p2)² = (E1 + E2)² - |p1 + p2|²
#[must_use]
#[inline]
pub fn invariant_mass_two_body(p1: &FourMomentum, p2: &FourMomentum) -> f64 {
    let total = *p1 + *p2;
    total.invariant_mass()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn particle_at_rest() {
        let p = FourMomentum::at_rest(938.272); // proton
        assert!((p.invariant_mass() - 938.272).abs() < 0.001);
        assert!((p.beta()).abs() < 1e-10);
        assert!((p.gamma() - 1.0).abs() < 1e-10);
        assert!((p.kinetic_energy()).abs() < 0.001);
    }

    #[test]
    fn massless_particle() {
        let photon = FourMomentum::new(1000.0, 0.0, 0.0, 1000.0);
        assert!(photon.invariant_mass() < 1e-10);
        assert!((photon.beta() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn lorentz_gamma_at_rest() {
        assert!((lorentz_gamma(0.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn lorentz_gamma_half_c() {
        let g = lorentz_gamma(0.5);
        let expected = 1.0 / libm::sqrt(0.75);
        assert!((g - expected).abs() < 1e-10);
    }

    #[test]
    fn lorentz_gamma_at_c() {
        assert!(lorentz_gamma(1.0).is_infinite());
    }

    #[test]
    fn gamma_beta_roundtrip() {
        let beta = 0.8;
        let gamma = lorentz_gamma(beta);
        let beta_back = gamma_to_beta(gamma);
        assert!((beta - beta_back).abs() < 1e-10);
    }

    #[test]
    fn velocity_addition_subluminal() {
        // 0.5c + 0.5c = 0.8c (not 1.0c)
        let result = velocity_addition(0.5, 0.5);
        assert!((result - 0.8).abs() < 1e-10);
    }

    #[test]
    fn velocity_addition_with_c() {
        // anything + c = c
        let result = velocity_addition(0.9, 1.0);
        assert!((result - 1.0).abs() < 1e-10);
    }

    #[test]
    fn energy_momentum_relation() {
        // E² = (mc²)² + (pc)²
        let mass = 938.272; // proton MeV
        let momentum = 500.0; // MeV/c
        let energy = relativistic_energy(mass, momentum);
        let e_check = libm::sqrt(mass * mass + momentum * momentum);
        assert!((energy - e_check).abs() < 1e-6);
    }

    #[test]
    fn momentum_from_kinetic_energy() {
        let mass = 0.511; // electron MeV
        let ke = 1.0; // 1 MeV kinetic energy
        let p = relativistic_momentum(mass, ke);
        // E = T + m = 1.511, p = sqrt(E² - m²) = sqrt(1.511² - 0.511²)
        let expected = libm::sqrt(1.511 * 1.511 - 0.511 * 0.511);
        assert!((p - expected).abs() < 0.001, "p={p}, expected={expected}");
    }

    #[test]
    fn de_broglie_proton() {
        // 1 GeV/c proton: λ = 2π × 197.3 / 1000 ≈ 1.24 fm
        let lambda = de_broglie_wavelength_fm(1000.0);
        assert!(lambda > 1.0 && lambda < 2.0, "λ={lambda} fm");
    }

    #[test]
    fn invariant_mass_two_body_at_rest() {
        let p1 = FourMomentum::at_rest(938.272);
        let p2 = FourMomentum::at_rest(938.272);
        let m_inv = invariant_mass_two_body(&p1, &p2);
        // Two protons at rest: M = 2 × 938.272
        assert!((m_inv - 2.0 * 938.272).abs() < 0.01);
    }

    #[test]
    fn four_momentum_add() {
        let p1 = FourMomentum::new(100.0, 10.0, 0.0, 50.0);
        let p2 = FourMomentum::new(200.0, -10.0, 5.0, 30.0);
        let total = p1 + p2;
        assert!((total.energy - 300.0).abs() < 1e-10);
        assert!((total.px).abs() < 1e-10);
        assert!((total.py - 5.0).abs() < 1e-10);
        assert!((total.pz - 80.0).abs() < 1e-10);
    }

    #[test]
    fn serde_roundtrip_four_momentum() {
        let p = FourMomentum::new(938.272, 100.0, 200.0, 300.0);
        let json = serde_json::to_string(&p).unwrap();
        let back: FourMomentum = serde_json::from_str(&json).unwrap();
        assert_eq!(p, back);
    }

    #[test]
    fn from_mass_and_momentum() {
        let p = FourMomentum::from_mass_and_momentum(938.272, 500.0);
        assert!(p.energy > 938.272); // E > m for nonzero p
        assert!((p.pz - 500.0).abs() < 1e-10);
        assert!((p.invariant_mass() - 938.272).abs() < 0.01);
    }

    #[test]
    fn massless_gamma_infinite() {
        let photon = FourMomentum::new(100.0, 0.0, 0.0, 100.0);
        assert!(photon.gamma().is_infinite());
    }

    #[test]
    fn gamma_to_beta_at_rest() {
        assert!((gamma_to_beta(1.0)).abs() < 1e-10);
    }

    #[test]
    fn velocity_conversions() {
        let v = 1.5e8; // m/s
        let beta = velocity_to_beta(v);
        let v_back = beta_to_velocity(beta);
        assert!((v - v_back).abs() < 1.0);
    }

    #[test]
    fn de_broglie_zero_momentum() {
        assert!(de_broglie_wavelength_fm(0.0).is_infinite());
    }

    #[test]
    fn invariant_mass_ultra_relativistic_electron() {
        // 10 GeV/c electron: E² − p² would lose ~6 digits; (E − p)(E + p) keeps them.
        let p = FourMomentum::from_mass_and_momentum(0.510_998_950_69, 1e4);
        assert!((p.invariant_mass() - 0.510_998_950_69).abs() / 0.511 < 1e-6);
        assert!((p.kinetic_energy() - (p.energy - 0.510_998_950_69)).abs() < 1e-7);
    }

    #[test]
    fn gamma_symmetric_and_edge_cases() {
        assert!((lorentz_gamma(-0.6) - 1.25).abs() < 1e-12);
        assert!(lorentz_gamma(-1.5).is_infinite());
        assert!(velocity_addition(1.0, -1.0).is_nan());
        let beta = 1.0 - 1e-12;
        assert!((gamma_to_beta(lorentz_gamma(beta)) - beta).abs() < 1e-15);
    }

    #[test]
    fn zero_energy_beta() {
        let p = FourMomentum::new(0.0, 0.0, 0.0, 0.0);
        assert!((p.beta()).abs() < 1e-10);
    }
}
