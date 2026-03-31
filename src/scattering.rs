//! Scattering theory: Rutherford, Mott, Born, Compton, and pair production cross-sections.
//!
//! Provides differential and total cross-section calculations for
//! charged particle scattering off nuclei, partial-wave analysis,
//! Klein-Nishina Compton scattering, and Bethe-Heitler pair production.

use crate::constants::{
    AMU_MEV, CLASSICAL_ELECTRON_RADIUS_FM, ELECTRON_MASS_MEV, FINE_STRUCTURE, PROTON_MASS_MEV,
};

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

// ─── Partial-Wave Analysis ───────────────────────────────────────────────────

/// Computes the Legendre polynomial P_l(x) using the standard recurrence.
///
/// P_0(x) = 1, P_1(x) = x,
/// (l+1) P_{l+1}(x) = (2l+1) x P_l(x) - l P_{l-1}(x)
///
/// `x` should be in [-1, 1] (typically cos θ).
#[must_use]
#[inline]
pub fn legendre_polynomial(l: u32, x: f64) -> f64 {
    match l {
        0 => 1.0,
        1 => x,
        _ => {
            let mut p_prev = 1.0_f64;
            let mut p_curr = x;
            for n in 1..l {
                let n_f = n as f64;
                let p_next = ((2.0 * n_f + 1.0) * x * p_curr - n_f * p_prev) / (n_f + 1.0);
                p_prev = p_curr;
                p_curr = p_next;
            }
            p_curr
        }
    }
}

/// Total cross-section from partial-wave phase shifts in fm².
///
/// σ_total = (4π/k²) Σ_l (2l+1) sin²(δ_l)
///
/// Parameters:
/// - `k_inv_fm`: wavenumber in fm⁻¹
/// - `phase_shifts`: slice where `phase_shifts[l]` is δ_l in radians
#[must_use]
pub fn partial_wave_cross_section(k_inv_fm: f64, phase_shifts: &[f64]) -> f64 {
    if k_inv_fm <= 0.0 {
        return 0.0;
    }
    let k2 = k_inv_fm * k_inv_fm;
    let mut sum = 0.0_f64;
    for (l, &delta) in phase_shifts.iter().enumerate() {
        let sin_d = libm::sin(delta);
        sum += (2.0 * l as f64 + 1.0) * sin_d * sin_d;
    }
    4.0 * core::f64::consts::PI * sum / k2
}

/// Differential cross-section |f(θ)|² from partial-wave phase shifts in fm²/sr.
///
/// f(θ) = (1/k) Σ_l (2l+1) exp(iδ_l) sin(δ_l) P_l(cos θ)
///
/// Parameters:
/// - `k_inv_fm`: wavenumber in fm⁻¹
/// - `phase_shifts`: slice where `phase_shifts[l]` is δ_l in radians
/// - `theta_rad`: scattering angle in radians
#[must_use]
pub fn partial_wave_differential(k_inv_fm: f64, phase_shifts: &[f64], theta_rad: f64) -> f64 {
    if k_inv_fm <= 0.0 {
        return 0.0;
    }
    let cos_theta = libm::cos(theta_rad);
    let mut re_f = 0.0_f64;
    let mut im_f = 0.0_f64;
    for (l, &delta) in phase_shifts.iter().enumerate() {
        let sin_d = libm::sin(delta);
        let cos_d = libm::cos(delta);
        let pl = legendre_polynomial(l as u32, cos_theta);
        let weight = (2.0 * l as f64 + 1.0) * pl;
        // exp(iδ) sin(δ) = (cos δ + i sin δ) sin δ = sin δ cos δ + i sin²δ
        re_f += weight * sin_d * cos_d;
        im_f += weight * sin_d * sin_d;
    }
    (re_f * re_f + im_f * im_f) / (k_inv_fm * k_inv_fm)
}

// ─── Born Approximation for Screened Coulomb ─────────────────────────────────

/// Thomas-Fermi screening length in femtometers.
///
/// a_TF = 0.8853 × a₀ / Z^(1/3)
///
/// where a₀ = 52917.72 fm (Bohr radius in fm).
///
/// Reference: Thomas-Fermi model, see e.g. Ziegler, Biersack & Littmark (1985).
#[must_use]
#[inline]
pub fn thomas_fermi_screening_fm(z: u32) -> f64 {
    if z == 0 {
        return 0.0;
    }
    // Bohr radius in fm: 5.29177210903e-11 m * 1e15 fm/m = 52917.72 fm
    let a0_fm = 52_917.72;
    0.8853 * a0_fm / libm::cbrt(z as f64)
}

/// Born approximation differential cross-section for screened Coulomb (Yukawa)
/// potential in fm²/sr.
///
/// V(r) = (Z₁ Z₂ α ħc / r) exp(-r/a)
///
/// dσ/dΩ = [2μ Z₁ Z₂ α ħc / (ħ² (q² + 1/a²))]²
///       = [Z₁ Z₂ α ħc / (2E_cm)]² × 1 / [sin²(θ/2) + (ħ/(2ka))²]²
///       (after simplification with ħ²k² = 2μE)
///
/// Parameters:
/// - `z_proj`, `z_target`: atomic numbers
/// - `energy_cm_mev`: center-of-mass kinetic energy in MeV
/// - `theta_rad`: scattering angle in radians
/// - `screening_fm`: screening length a in fm
#[must_use]
pub fn born_screened_coulomb(
    z_proj: u32,
    z_target: u32,
    energy_cm_mev: f64,
    theta_rad: f64,
    screening_fm: f64,
) -> f64 {
    if energy_cm_mev <= 0.0 || screening_fm <= 0.0 {
        return 0.0;
    }

    let hbar_c = 197.326_980_4; // MeV·fm

    // Reduced mass approximation: A ~ 2Z for Z > 1, proton mass for Z = 1.
    let mass_proj = if z_proj <= 1 {
        PROTON_MASS_MEV
    } else {
        z_proj as f64 * 2.0 * AMU_MEV
    };
    let mass_target = if z_target <= 1 {
        PROTON_MASS_MEV
    } else {
        z_target as f64 * 2.0 * AMU_MEV
    };
    let mu = mass_proj * mass_target / (mass_proj + mass_target);

    let sin_half = libm::sin(theta_rad / 2.0);
    let sin2_half = sin_half * sin_half;

    // k² = 2μE/(ħc)²
    let k2 = 2.0 * mu * energy_cm_mev / (hbar_c * hbar_c);
    // q² = 4k² sin²(θ/2)
    let q2 = 4.0 * k2 * sin2_half;
    // Screening: q² + 1/a²
    let inv_a2 = 1.0 / (screening_fm * screening_fm);
    let denom = q2 + inv_a2;
    if denom < 1e-30 {
        return f64::INFINITY;
    }

    // dσ/dΩ = (2μ Z1 Z2 α / ((ħc)(q²+1/a²)))²
    let numer = 2.0 * mu * z_proj as f64 * z_target as f64 * FINE_STRUCTURE / hbar_c;
    (numer * numer) / (denom * denom)
}

// ─── Electron-Atom Elastic Scattering (Mott with Form Factor) ────────────────

/// Mott differential cross-section for relativistic electrons on a point nucleus.
///
/// dσ/dΩ = (Z α ħc / (2E))² × cos²(θ/2) / sin⁴(θ/2)
///
/// Parameters:
/// - `z_target`: target atomic number
/// - `electron_energy_mev`: electron kinetic energy in MeV (relativistic)
/// - `theta_rad`: scattering angle in radians
///
/// Returns dσ/dΩ in fm²/sr.
#[must_use]
#[inline]
pub fn mott_electron_differential(z_target: u32, electron_energy_mev: f64, theta_rad: f64) -> f64 {
    if electron_energy_mev <= 0.0 {
        return 0.0;
    }
    let hbar_c = 197.326_980_4;
    let sin_half = libm::sin(theta_rad / 2.0);
    let cos_half = libm::cos(theta_rad / 2.0);
    let sin4 = sin_half * sin_half * sin_half * sin_half;
    if sin4 < 1e-30 {
        return f64::INFINITY;
    }
    let a_param = z_target as f64 * FINE_STRUCTURE * hbar_c / (2.0 * electron_energy_mev);
    a_param * a_param * cos_half * cos_half / sin4
}

/// Nuclear form factor for a uniform charge distribution (sphere of radius R).
///
/// F(q) = 3 [sin(qR) - qR cos(qR)] / (qR)³
///
/// Returns 1.0 for qR < 1e-6 to avoid division by zero.
///
/// Parameters:
/// - `q_inv_fm`: momentum transfer in fm⁻¹
/// - `radius_fm`: nuclear charge radius in fm
#[must_use]
#[inline]
pub fn nuclear_form_factor_uniform(q_inv_fm: f64, radius_fm: f64) -> f64 {
    let qr = q_inv_fm * radius_fm;
    if libm::fabs(qr) < 1e-6 {
        return 1.0;
    }
    let sin_qr = libm::sin(qr);
    let cos_qr = libm::cos(qr);
    let qr3 = qr * qr * qr;
    3.0 * (sin_qr - qr * cos_qr) / qr3
}

/// Mott electron scattering cross-section with nuclear form factor.
///
/// dσ/dΩ = dσ/dΩ_Mott × |F(q)|²
///
/// Uses nuclear radius R = 1.2 × A^(1/3) fm and momentum transfer
/// q = 2E sin(θ/2) / (ħc) for relativistic electrons.
///
/// Parameters:
/// - `z_target`: target atomic number
/// - `a_target`: target mass number
/// - `electron_energy_mev`: electron energy in MeV
/// - `theta_rad`: scattering angle in radians
///
/// Returns dσ/dΩ in fm²/sr.
#[must_use]
pub fn mott_electron_with_form_factor(
    z_target: u32,
    a_target: u32,
    electron_energy_mev: f64,
    theta_rad: f64,
) -> f64 {
    let hbar_c = 197.326_980_4;
    let mott = mott_electron_differential(z_target, electron_energy_mev, theta_rad);

    // Momentum transfer for relativistic electrons: q = 2E sin(θ/2) / ħc
    let sin_half = libm::sin(theta_rad / 2.0);
    let q = 2.0 * electron_energy_mev * sin_half / hbar_c;

    // Nuclear radius: R = 1.2 * A^(1/3) fm
    let radius = 1.2 * libm::cbrt(a_target as f64);

    let ff = nuclear_form_factor_uniform(q, radius);
    mott * ff * ff
}

// ─── Compton Scattering (Klein-Nishina) ──────────────────────────────────────

/// Energy ratio ω'/ω for Compton scattering.
///
/// ω'/ω = 1 / [1 + (E_γ / m_e c²)(1 - cos θ)]
///
/// Parameters:
/// - `photon_energy_mev`: incident photon energy in MeV
/// - `theta_rad`: scattering angle in radians
#[must_use]
#[inline]
pub fn compton_energy_ratio(photon_energy_mev: f64, theta_rad: f64) -> f64 {
    let gamma = photon_energy_mev / ELECTRON_MASS_MEV;
    let cos_theta = libm::cos(theta_rad);
    1.0 / (1.0 + gamma * (1.0 - cos_theta))
}

/// Klein-Nishina differential cross-section for Compton scattering in fm²/sr.
///
/// dσ/dΩ = (r_e²/2) (ω'/ω)² [ω/ω' + ω'/ω - sin²θ]
///
/// Parameters:
/// - `photon_energy_mev`: incident photon energy in MeV
/// - `theta_rad`: scattering angle in radians
#[must_use]
pub fn klein_nishina_differential(photon_energy_mev: f64, theta_rad: f64) -> f64 {
    if photon_energy_mev <= 0.0 {
        return 0.0;
    }
    let ratio = compton_energy_ratio(photon_energy_mev, theta_rad);
    let inv_ratio = 1.0 / ratio;
    let sin_theta = libm::sin(theta_rad);
    let re2_half = CLASSICAL_ELECTRON_RADIUS_FM * CLASSICAL_ELECTRON_RADIUS_FM / 2.0;
    re2_half * ratio * ratio * (inv_ratio + ratio - sin_theta * sin_theta)
}

/// Total Klein-Nishina cross-section for Compton scattering in fm².
///
/// σ_KN = 2π r_e² { [(1+γ)/γ³] [2γ(1+γ)/(1+2γ) - ln(1+2γ)]
///         + ln(1+2γ)/(2γ) - (1+3γ)/(1+2γ)² }
///
/// where γ = E_γ / (m_e c²).
///
/// In the low-energy limit (γ → 0), returns the Thomson cross-section
/// σ_T = (8/3) π r_e².
#[must_use]
pub fn klein_nishina_total(photon_energy_mev: f64) -> f64 {
    if photon_energy_mev <= 0.0 {
        return 0.0;
    }
    let re2 = CLASSICAL_ELECTRON_RADIUS_FM * CLASSICAL_ELECTRON_RADIUS_FM;
    let gamma = photon_energy_mev / ELECTRON_MASS_MEV;

    // For very low photon energies, use Thomson limit to avoid numerical issues
    if gamma < 1e-6 {
        return 8.0 * core::f64::consts::PI * re2 / 3.0;
    }

    let g2 = gamma * gamma;
    let g3 = g2 * gamma;
    let one_plus_2g = 1.0 + 2.0 * gamma;
    let ln_term = libm::log(one_plus_2g);

    let term1 = (1.0 + gamma) / g3 * (2.0 * gamma * (1.0 + gamma) / one_plus_2g - ln_term);
    let term2 = ln_term / (2.0 * gamma);
    let term3 = (1.0 + 3.0 * gamma) / (one_plus_2g * one_plus_2g);

    2.0 * core::f64::consts::PI * re2 * (term1 + term2 - term3)
}

// ─── Pair Production (Bethe-Heitler) ─────────────────────────────────────────

/// Bethe-Heitler pair production cross-section (high-energy asymptotic limit).
///
/// σ_pair ≈ α r_e² Z² [28/9 ln(2E/(m_e c²)) - 218/27]
///
/// Valid for E_γ >> 2 m_e c² = 1.022 MeV. Returns 0.0 below threshold.
///
/// Parameters:
/// - `z_target`: target atomic number
/// - `photon_energy_mev`: photon energy in MeV
///
/// Returns cross-section in fm² per atom.
#[must_use]
pub fn pair_production_cross_section(z_target: u32, photon_energy_mev: f64) -> f64 {
    let threshold = 2.0 * ELECTRON_MASS_MEV; // 1.021998 MeV
    if photon_energy_mev <= threshold {
        return 0.0;
    }
    let re2 = CLASSICAL_ELECTRON_RADIUS_FM * CLASSICAL_ELECTRON_RADIUS_FM;
    let z2 = (z_target as f64) * (z_target as f64);
    let ratio = 2.0 * photon_energy_mev / ELECTRON_MASS_MEV;
    let ln_ratio = libm::log(ratio);

    FINE_STRUCTURE * re2 * z2 * (28.0 / 9.0 * ln_ratio - 218.0 / 27.0)
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

    // ─── Legendre Polynomial Tests ───────────────────────────────────────────

    #[test]
    fn legendre_p0_is_one() {
        assert!((legendre_polynomial(0, 0.5) - 1.0).abs() < 1e-12);
        assert!((legendre_polynomial(0, -0.3) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn legendre_p1_is_x() {
        let x = 0.7;
        assert!((legendre_polynomial(1, x) - x).abs() < 1e-12);
    }

    #[test]
    fn legendre_p2_formula() {
        // P_2(x) = (3x² - 1) / 2
        let x = 0.6;
        let expected = (3.0 * x * x - 1.0) / 2.0;
        assert!((legendre_polynomial(2, x) - expected).abs() < 1e-12);
    }

    #[test]
    fn legendre_p3_formula() {
        // P_3(x) = (5x³ - 3x) / 2
        let x = 0.4;
        let expected = (5.0 * x * x * x - 3.0 * x) / 2.0;
        assert!(
            (legendre_polynomial(3, x) - expected).abs() < 1e-12,
            "P_3({x}) = {}, expected {expected}",
            legendre_polynomial(3, x)
        );
    }

    #[test]
    fn legendre_at_one() {
        // P_l(1) = 1 for all l
        for l in 0..10 {
            assert!(
                (legendre_polynomial(l, 1.0) - 1.0).abs() < 1e-10,
                "P_{l}(1) should be 1"
            );
        }
    }

    // ─── Partial-Wave Analysis Tests ─────────────────────────────────────────

    #[test]
    fn partial_wave_s_wave_only() {
        // Only l=0 with δ₀ = π/4: σ = 4π/k² sin²(π/4) = 4π/k² * 0.5 = 2π/k²
        let k = 0.5;
        let sigma = partial_wave_cross_section(k, &[core::f64::consts::FRAC_PI_4]);
        let expected = 2.0 * core::f64::consts::PI / (k * k);
        assert!(
            (sigma - expected).abs() / expected < 1e-10,
            "σ={sigma}, expected={expected}"
        );
    }

    #[test]
    fn partial_wave_zero_phase_shifts() {
        // All δ_l = 0 -> σ = 0
        let sigma = partial_wave_cross_section(1.0, &[0.0, 0.0, 0.0]);
        assert!(sigma.abs() < 1e-15);
    }

    #[test]
    fn partial_wave_unitarity_limit() {
        // Unitarity limit for l=0: max σ from s-wave is 4π/k² (when δ₀ = π/2)
        let k = 1.0;
        let sigma = partial_wave_cross_section(k, &[core::f64::consts::FRAC_PI_2]);
        let max_s_wave = 4.0 * core::f64::consts::PI / (k * k);
        assert!(
            (sigma - max_s_wave).abs() / max_s_wave < 1e-10,
            "Unitarity limit: σ={sigma}, max={max_s_wave}"
        );
    }

    #[test]
    fn partial_wave_differential_isotropic_s_wave() {
        // Pure s-wave (l=0) scattering should be isotropic
        let k = 1.0;
        let phase = [core::f64::consts::FRAC_PI_4];
        let ds_0 = partial_wave_differential(k, &phase, 0.1);
        let ds_90 = partial_wave_differential(k, &phase, core::f64::consts::FRAC_PI_2);
        let ds_180 = partial_wave_differential(k, &phase, core::f64::consts::PI);
        // Should all be equal for isotropic
        assert!(
            (ds_0 - ds_90).abs() / ds_0 < 1e-10,
            "s-wave not isotropic: {ds_0} vs {ds_90}"
        );
        assert!(
            (ds_0 - ds_180).abs() / ds_0 < 1e-10,
            "s-wave not isotropic: {ds_0} vs {ds_180}"
        );
    }

    #[test]
    fn partial_wave_differential_positive() {
        let k = 0.5;
        let phase = [0.3, 0.1, 0.05];
        let ds = partial_wave_differential(k, &phase, 1.0);
        assert!(ds > 0.0 && ds.is_finite(), "dσ/dΩ={ds}");
    }

    // ─── Born Screened Coulomb Tests ─────────────────────────────────────────

    #[test]
    fn born_screened_approaches_rutherford_for_large_screening() {
        // Large screening length -> approaches Rutherford
        let theta = 1.0;
        let energy = 5.0;
        let ruth = rutherford_differential(1, 79, energy, theta);
        let born = born_screened_coulomb(1, 79, energy, theta, 1e6);
        let ratio = born / ruth;
        // Should be close to 1 for very large screening length
        assert!(
            (ratio - 1.0).abs() < 0.05,
            "Born/Ruth ratio={ratio}, should approach 1"
        );
    }

    #[test]
    fn born_screened_less_than_rutherford() {
        // Screening always reduces the cross-section vs Rutherford at forward angles
        let theta = 0.1;
        let energy = 5.0;
        let ruth = rutherford_differential(1, 79, energy, theta);
        let born = born_screened_coulomb(1, 79, energy, theta, 10000.0);
        assert!(
            born < ruth,
            "Screened ({born}) should be < Rutherford ({ruth}) at forward angles"
        );
    }

    #[test]
    fn born_screened_finite_at_zero_angle() {
        // Unlike Rutherford, screened Coulomb is finite at θ=0
        let ds = born_screened_coulomb(2, 79, 5.0, 0.0, 10000.0);
        assert!(ds.is_finite(), "Screened should be finite at θ=0: {ds}");
    }

    // ─── Thomas-Fermi Screening Tests ────────────────────────────────────────

    #[test]
    fn thomas_fermi_decreases_with_z() {
        // Screening length decreases with Z (a ∝ Z^(-1/3))
        let a_low = thomas_fermi_screening_fm(10);
        let a_high = thomas_fermi_screening_fm(80);
        assert!(
            a_high < a_low,
            "a(Z=80)={a_high} should be < a(Z=10)={a_low}"
        );
    }

    #[test]
    fn thomas_fermi_hydrogen() {
        // Z=1: a_TF = 0.8853 * a0 ≈ 0.8853 * 52917.72 ≈ 46,847 fm
        let a = thomas_fermi_screening_fm(1);
        let expected = 0.8853 * 52_917.72;
        assert!(
            (a - expected).abs() / expected < 1e-6,
            "a(Z=1)={a}, expected≈{expected}"
        );
    }

    #[test]
    fn thomas_fermi_zero_z_returns_zero() {
        assert!((thomas_fermi_screening_fm(0)).abs() < 1e-15);
    }

    // ─── Mott Electron Scattering Tests ──────────────────────────────────────

    #[test]
    fn mott_electron_diverges_forward() {
        let ds = mott_electron_differential(79, 200.0, 0.001);
        assert!(ds > 1e10, "Forward Mott should be very large: {ds}");
    }

    #[test]
    fn mott_electron_finite_at_90() {
        let ds = mott_electron_differential(79, 200.0, core::f64::consts::FRAC_PI_2);
        assert!(ds > 0.0 && ds.is_finite(), "90° Mott = {ds}");
    }

    #[test]
    fn mott_electron_scales_with_z_squared() {
        let ds1 = mott_electron_differential(40, 200.0, 1.0);
        let ds2 = mott_electron_differential(80, 200.0, 1.0);
        let ratio = ds2 / ds1;
        assert!((ratio - 4.0).abs() < 0.01, "Z² scaling: ratio={ratio}");
    }

    // ─── Nuclear Form Factor Tests ───────────────────────────────────────────

    #[test]
    fn form_factor_unity_at_zero_q() {
        assert!((nuclear_form_factor_uniform(0.0, 5.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn form_factor_less_than_one_nonzero_q() {
        let ff = nuclear_form_factor_uniform(0.5, 5.0);
        assert!(
            libm::fabs(ff) < 1.0,
            "F(q) should be < 1 for finite q: {ff}"
        );
    }

    #[test]
    fn form_factor_has_zeros() {
        // The uniform sphere form factor has zeros at tan(qR) = qR
        // First zero is at qR ≈ 4.493
        let r = 5.0;
        let q = 4.493 / r;
        let ff = nuclear_form_factor_uniform(q, r);
        assert!(
            libm::fabs(ff) < 0.01,
            "F(q) near first zero should be ~0: {ff}"
        );
    }

    // ─── Mott with Form Factor Tests ─────────────────────────────────────────

    #[test]
    fn mott_with_ff_less_than_point() {
        // Form factor always reduces cross-section at finite angles
        let point = mott_electron_differential(82, 200.0, 0.5);
        let with_ff = mott_electron_with_form_factor(82, 208, 200.0, 0.5);
        assert!(
            with_ff <= point,
            "With FF ({with_ff}) should be <= point ({point})"
        );
    }

    #[test]
    fn mott_with_ff_approaches_point_at_small_angle() {
        // At very small angles (small q), F(q) ≈ 1
        let point = mott_electron_differential(6, 10.0, 0.01);
        let with_ff = mott_electron_with_form_factor(6, 12, 10.0, 0.01);
        let ratio = with_ff / point;
        assert!((ratio - 1.0).abs() < 0.01, "Small angle: ratio={ratio}");
    }

    #[test]
    fn mott_with_ff_positive() {
        let ds = mott_electron_with_form_factor(79, 197, 500.0, 0.3);
        assert!(ds >= 0.0 && ds.is_finite(), "dσ/dΩ={ds}");
    }

    // ─── Compton Scattering Tests ────────────────────────────────────────────

    #[test]
    fn compton_ratio_unity_at_zero_angle() {
        // θ=0: no energy loss, ω'/ω = 1
        let ratio = compton_energy_ratio(1.0, 0.0);
        assert!((ratio - 1.0).abs() < 1e-12);
    }

    #[test]
    fn compton_ratio_backscatter() {
        // θ=π: maximum energy loss, ω'/ω = 1/(1+2γ)
        let energy = 0.511; // ~1 m_e c²
        let gamma = energy / ELECTRON_MASS_MEV;
        let ratio = compton_energy_ratio(energy, core::f64::consts::PI);
        let expected = 1.0 / (1.0 + 2.0 * gamma);
        assert!(
            (ratio - expected).abs() < 1e-10,
            "ratio={ratio}, expected={expected}"
        );
    }

    #[test]
    fn compton_ratio_decreases_with_energy() {
        let r_low = compton_energy_ratio(0.1, 1.0);
        let r_high = compton_energy_ratio(10.0, 1.0);
        assert!(r_high < r_low, "Higher energy should lose more");
    }

    #[test]
    fn klein_nishina_differential_positive() {
        let ds = klein_nishina_differential(0.662, core::f64::consts::FRAC_PI_2);
        assert!(ds > 0.0 && ds.is_finite(), "dσ/dΩ={ds}");
    }

    #[test]
    fn klein_nishina_low_energy_approaches_thomson() {
        // At low energy, KN → Thomson: dσ/dΩ = r_e²/2 * (1 + cos²θ)
        let energy = 1e-4; // very low
        let theta = core::f64::consts::FRAC_PI_2;
        let ds = klein_nishina_differential(energy, theta);
        let re2 = CLASSICAL_ELECTRON_RADIUS_FM * CLASSICAL_ELECTRON_RADIUS_FM;
        let cos_t = libm::cos(theta);
        let thomson = re2 / 2.0 * (1.0 + cos_t * cos_t);
        assert!(
            (ds - thomson).abs() / thomson < 1e-3,
            "KN={ds}, Thomson={thomson}"
        );
    }

    #[test]
    fn klein_nishina_forward_backward_symmetry_low_energy() {
        // At low energy, forward and backward should be similar
        let energy = 1e-4;
        let ds_fwd = klein_nishina_differential(energy, 0.3);
        let ds_bwd = klein_nishina_differential(energy, core::f64::consts::PI - 0.3);
        let ratio = ds_fwd / ds_bwd;
        assert!((ratio - 1.0).abs() < 0.01, "Low-energy asymmetry: {ratio}");
    }

    #[test]
    fn klein_nishina_total_thomson_limit() {
        // σ_T = (8/3) π r_e² ≈ 66.52 fm²
        let re2 = CLASSICAL_ELECTRON_RADIUS_FM * CLASSICAL_ELECTRON_RADIUS_FM;
        let thomson = 8.0 * core::f64::consts::PI * re2 / 3.0;
        let sigma = klein_nishina_total(1e-4);
        assert!(
            (sigma - thomson).abs() / thomson < 1e-3,
            "σ_KN={sigma}, σ_T={thomson}"
        );
    }

    #[test]
    fn klein_nishina_total_decreases_with_energy() {
        let s1 = klein_nishina_total(0.1);
        let s2 = klein_nishina_total(1.0);
        let s3 = klein_nishina_total(10.0);
        assert!(s1 > s2 && s2 > s3, "σ should decrease: {s1}, {s2}, {s3}");
    }

    #[test]
    fn klein_nishina_total_positive() {
        let sigma = klein_nishina_total(0.662);
        assert!(sigma > 0.0 && sigma.is_finite(), "σ={sigma}");
    }

    // ─── Pair Production Tests ───────────────────────────────────────────────

    #[test]
    fn pair_production_zero_below_threshold() {
        // Below 1.022 MeV, no pair production
        assert!((pair_production_cross_section(82, 0.5)).abs() < 1e-15);
        assert!((pair_production_cross_section(82, 1.02)).abs() < 1e-15);
    }

    #[test]
    fn pair_production_positive_above_threshold() {
        let sigma = pair_production_cross_section(82, 10.0);
        assert!(sigma > 0.0 && sigma.is_finite(), "σ={sigma}");
    }

    #[test]
    fn pair_production_scales_with_z_squared() {
        let s1 = pair_production_cross_section(40, 100.0);
        let s2 = pair_production_cross_section(80, 100.0);
        let ratio = s2 / s1;
        assert!((ratio - 4.0).abs() < 0.01, "Z² scaling: ratio={ratio}");
    }

    #[test]
    fn pair_production_increases_with_energy() {
        let s10 = pair_production_cross_section(82, 10.0);
        let s100 = pair_production_cross_section(82, 100.0);
        let s1000 = pair_production_cross_section(82, 1000.0);
        assert!(
            s100 > s10 && s1000 > s100,
            "σ should increase with energy: {s10}, {s100}, {s1000}"
        );
    }
}
