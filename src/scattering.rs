//! Scattering theory: Rutherford and Mott cross-sections.
//!
//! Provides differential and total cross-section calculations for
//! charged particle scattering off nuclei.

use crate::constants::FINE_STRUCTURE;

/// Calculates the Rutherford differential cross-section dσ/dΩ in fm².
///
/// The Rutherford formula for Coulomb scattering of a point charge off
/// a point nucleus:
///
/// dσ/dΩ = (Z₁ Z₂ e² / 4E)² / sin⁴(θ/2)
///
/// Using natural units: dσ/dΩ = (Z₁ Z₂ α ħc / (4 T_cm))² / sin⁴(θ/2)
///
/// Parameters:
/// - `z_proj`: atomic number of projectile
/// - `z_target`: atomic number of target
/// - `energy_cm_mev`: center-of-mass kinetic energy in MeV
/// - `theta_rad`: scattering angle in radians
///
/// Returns differential cross-section in fm²/sr.
#[must_use]
#[inline]
pub fn rutherford_differential(
    z_proj: u32,
    z_target: u32,
    energy_cm_mev: f64,
    theta_rad: f64,
) -> f64 {
    if energy_cm_mev <= 0.0 {
        return 0.0;
    }

    let sin_half = libm::sin(theta_rad / 2.0);
    let sin4 = sin_half * sin_half * sin_half * sin_half;
    if sin4 < 1e-30 {
        return f64::INFINITY; // forward scattering divergence
    }

    // a = Z₁ Z₂ α ħc / (4 T_cm), where ħc = 197.3269804 MeV·fm
    let hbar_c = 197.326_980_4;
    let a_param = z_proj as f64 * z_target as f64 * FINE_STRUCTURE * hbar_c / (4.0 * energy_cm_mev);

    (a_param * a_param) / sin4
}

/// Calculates the Mott scattering correction factor for identical particles.
///
/// The Mott correction accounts for quantum mechanical exchange effects
/// when identical particles scatter (e.g., proton-proton). The cross-section
/// is modified to:
///
/// dσ/dΩ_Mott = dσ/dΩ_Ruth × [1 - β² sin²(θ/2)]
///
/// for spin-0 identical particles (approximation), where β = v/c.
///
/// Parameters:
/// - `beta`: v/c of the projectile in the CM frame
/// - `theta_rad`: scattering angle in radians
///
/// Returns the Mott correction factor (multiply by Rutherford result).
#[must_use]
#[inline]
pub fn mott_correction_factor(beta: f64, theta_rad: f64) -> f64 {
    let sin_half = libm::sin(theta_rad / 2.0);
    1.0 - beta * beta * sin_half * sin_half
}

/// Calculates the Rutherford cross-section integrated over a minimum angle.
///
/// σ(θ > θ_min) = π (a/2)² / tan²(θ_min/2)
///
/// where a = Z₁ Z₂ e² / (2 T_cm).
///
/// Returns total cross-section in fm² for scattering beyond θ_min.
#[must_use]
#[inline]
pub fn rutherford_total_above_angle(
    z_proj: u32,
    z_target: u32,
    energy_cm_mev: f64,
    theta_min_rad: f64,
) -> f64 {
    if energy_cm_mev <= 0.0 || theta_min_rad <= 0.0 {
        return f64::INFINITY;
    }

    let hbar_c = 197.326_980_4;
    let a_param = z_proj as f64 * z_target as f64 * FINE_STRUCTURE * hbar_c / (4.0 * energy_cm_mev);

    let tan_half = libm::tan(theta_min_rad / 2.0);
    let tan2 = tan_half * tan_half;
    if tan2 < 1e-30 {
        return f64::INFINITY;
    }

    core::f64::consts::PI * a_param * a_param / tan2
}

/// Calculates the distance of closest approach for Coulomb scattering (fm).
///
/// d = Z₁ Z₂ e² / (2 T_cm) = Z₁ Z₂ α ħc / (2 T_cm)
///
/// This is the classical turning point for a head-on collision (θ = π).
#[must_use]
#[inline]
pub fn distance_of_closest_approach(z_proj: u32, z_target: u32, energy_cm_mev: f64) -> f64 {
    if energy_cm_mev <= 0.0 {
        return f64::INFINITY;
    }
    let hbar_c = 197.326_980_4;
    z_proj as f64 * z_target as f64 * FINE_STRUCTURE * hbar_c / (2.0 * energy_cm_mev)
}

/// Calculates the Sommerfeld parameter η for Coulomb scattering.
///
/// η = Z₁ Z₂ e² / (ħv) = Z₁ Z₂ α / β
///
/// where β = v/c. Large η means classical regime (Rutherford valid).
#[must_use]
#[inline]
pub fn sommerfeld_parameter(z_proj: u32, z_target: u32, beta: f64) -> f64 {
    if beta <= 0.0 {
        return f64::INFINITY;
    }
    z_proj as f64 * z_target as f64 * FINE_STRUCTURE / beta
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rutherford_diverges_forward() {
        // θ → 0 should give very large (infinite) cross-section
        let ds = rutherford_differential(79, 2, 5.0, 0.001);
        assert!(ds > 1e10, "Forward Rutherford should be huge: {ds}");
    }

    #[test]
    fn rutherford_finite_at_90_degrees() {
        // 5 MeV alpha on gold at 90°
        let ds = rutherford_differential(2, 79, 5.0, core::f64::consts::FRAC_PI_2);
        assert!(ds > 0.0 && ds.is_finite(), "90° should be finite: {ds}");
    }

    #[test]
    fn rutherford_scales_with_z_squared() {
        // Double Z_target -> 4x cross-section
        let ds1 = rutherford_differential(2, 40, 10.0, 1.0);
        let ds2 = rutherford_differential(2, 80, 10.0, 1.0);
        let ratio = ds2 / ds1;
        assert!((ratio - 4.0).abs() < 0.01, "Z² scaling: ratio={ratio}");
    }

    #[test]
    fn rutherford_decreases_with_energy() {
        // Higher energy -> smaller cross-section (1/E²)
        let ds5 = rutherford_differential(2, 79, 5.0, 1.0);
        let ds10 = rutherford_differential(2, 79, 10.0, 1.0);
        let ratio = ds5 / ds10;
        assert!((ratio - 4.0).abs() < 0.1, "1/E² scaling: ratio={ratio}");
    }

    #[test]
    fn mott_correction_nonrelativistic() {
        // β ≈ 0 -> Mott factor ≈ 1
        let factor = mott_correction_factor(0.01, 1.0);
        assert!((factor - 1.0).abs() < 0.01);
    }

    #[test]
    fn mott_correction_relativistic() {
        // β = 0.9 at θ = π/2 -> factor = 1 - 0.81 * sin²(π/4)
        let factor = mott_correction_factor(0.9, core::f64::consts::FRAC_PI_2);
        let expected = 1.0 - 0.81 * 0.5;
        assert!((factor - expected).abs() < 0.01, "Mott={factor}");
    }

    #[test]
    fn closest_approach_gold_alpha() {
        // 5 MeV alpha on Au-197: d ≈ Z1*Z2*alpha*hbar_c/(2*T_cm)
        // ≈ 2*79*0.00730*197.3/(2*5) ≈ 2*79*1.44/10 ≈ 22.8 fm
        let dist = distance_of_closest_approach(2, 79, 5.0);
        assert!(dist > 20.0 && dist < 50.0, "d={dist} fm");
    }

    #[test]
    fn sommerfeld_large_means_classical() {
        // Slow alpha on heavy nucleus: η >> 1 -> classical
        let eta = sommerfeld_parameter(2, 79, 0.05);
        assert!(eta > 1.0, "η={eta}, should be >> 1 for classical");
    }

    #[test]
    fn rutherford_total_finite_above_angle() {
        let sigma = rutherford_total_above_angle(2, 79, 5.0, 0.1);
        assert!(sigma > 0.0 && sigma.is_finite(), "σ={sigma} fm²");
    }

    #[test]
    fn serde_roundtrip_not_needed() {
        // This module has no custom types needing serde roundtrip
        // (uses only f64 inputs/outputs and FourMomentum from relativity)
    }
}
