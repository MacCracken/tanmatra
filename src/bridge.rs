//! Cross-crate bridges — convert primitive values from other AGNOS science crates
//! into tanmatra atomic/nuclear physics parameters and vice versa.
//!
//! Always available — takes primitive values (f64), no science crate deps.

use alloc::vec::Vec;

// ── Bijli bridges (electromagnetism) ───────────────────────────────────────

/// Convert electron orbital energy (eV) to photon emission wavelength (nm).
///
/// λ = hc / E, where hc ≈ 1239.842 eV·nm.
#[must_use]
#[inline]
pub fn energy_to_wavelength_nm(energy_ev: f64) -> f64 {
    if energy_ev <= 0.0 {
        return 0.0;
    }
    1239.842 / energy_ev
}

/// Convert nuclear charge number (Z) to Coulomb field strength (V/m)
/// at distance r (m).
///
/// E = k × Z × e / r², where k = 8.9876e9, e = 1.602e-19 C.
#[must_use]
#[inline]
pub fn nuclear_charge_field(atomic_number: u32, distance_m: f64) -> f64 {
    const COULOMB_K: f64 = 8.987_551_792e9;
    const E_CHARGE: f64 = 1.602_176_634e-19;
    if distance_m <= 0.0 {
        return 0.0;
    }
    COULOMB_K * atomic_number as f64 * E_CHARGE / (distance_m * distance_m)
}

// ── Kimiya bridges (chemistry) ─────────────────────────────────────────────

/// Convert atomic number to number of valence electrons (main-group simplified).
#[must_use]
pub fn atomic_number_to_valence(atomic_number: u8) -> u8 {
    match atomic_number {
        1 => 1,
        2 => 2,
        3..=4 => atomic_number - 2,
        5..=10 => atomic_number - 2,
        11..=12 => atomic_number - 10,
        13..=18 => atomic_number - 10,
        _ => 2, // transition metals: conventional
    }
}

/// Convert atomic number and mass number to neutron count.
#[must_use]
#[inline]
pub fn neutron_count(atomic_number: u32, mass_number: u32) -> u32 {
    mass_number.saturating_sub(atomic_number)
}

/// Convert isotope mass (u) to nuclear binding energy deficit (MeV).
///
/// ΔE = (Z×m_p + N×m_n - M) × 931.494 MeV/u
/// `mass_u`: atomic mass in unified atomic mass units.
#[must_use]
pub fn mass_to_binding_deficit_mev(atomic_number: u32, mass_number: u32, mass_u: f64) -> f64 {
    const M_PROTON: f64 = 1.007_276_47;
    const M_NEUTRON: f64 = 1.008_664_92;
    const M_ELECTRON: f64 = 0.000_548_58;
    const MEV_PER_U: f64 = 931.494;

    let z = atomic_number as f64;
    let n = mass_number.saturating_sub(atomic_number) as f64;
    let constituents = z * (M_PROTON + M_ELECTRON) + n * M_NEUTRON;
    (constituents - mass_u) * MEV_PER_U
}

/// Convert half-life (seconds) to decay probability per unit time (1/s).
///
/// λ = ln(2) / t½
#[must_use]
#[inline]
pub fn half_life_to_decay_constant(half_life_s: f64) -> f64 {
    if half_life_s <= 0.0 {
        return 0.0;
    }
    core::f64::consts::LN_2 / half_life_s
}

// ── Prakash bridges (optics) ───────────────────────────────────────────────

/// Convert energy level transitions to spectral line wavelengths (nm).
///
/// Given a list of transition energies (eV), returns corresponding wavelengths.
#[must_use]
pub fn transitions_to_wavelengths(transition_energies_ev: &[f64]) -> Vec<f64> {
    transition_energies_ev
        .iter()
        .filter(|&&e| e > 0.0)
        .map(|&e| 1239.842 / e)
        .collect()
}

/// Convert nuclear spin quantum number to hyperfine splitting factor.
///
/// The hyperfine splitting scales with the nuclear magnetic moment,
/// which is roughly proportional to I (nuclear spin).
/// Returns a dimensionless scaling factor.
#[must_use]
#[inline]
pub fn nuclear_spin_to_hyperfine_scale(spin_i: f64) -> f64 {
    if spin_i <= 0.0 {
        return 0.0;
    }
    // Hyperfine splitting ∝ (2I+1) states
    2.0 * spin_i + 1.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn energy_to_wavelength_lyman_alpha() {
        // Lyman-α: ~10.2 eV → ~121.6 nm
        let nm = energy_to_wavelength_nm(10.2);
        assert!((nm - 121.6).abs() < 1.0);
    }

    #[test]
    fn energy_to_wavelength_zero() {
        assert!(energy_to_wavelength_nm(0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn nuclear_field_hydrogen() {
        // Z=1 at 1 Bohr radius → very strong field
        let e = nuclear_charge_field(1, 5.29e-11);
        assert!(e > 1e11, "got {e}");
    }

    #[test]
    fn nuclear_field_zero_distance() {
        assert!(nuclear_charge_field(1, 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn valence_carbon() {
        assert_eq!(atomic_number_to_valence(6), 4);
    }

    #[test]
    fn valence_chlorine() {
        assert_eq!(atomic_number_to_valence(17), 7);
    }

    #[test]
    fn neutron_count_carbon12() {
        assert_eq!(neutron_count(6, 12), 6);
    }

    #[test]
    fn neutron_count_underflow() {
        assert_eq!(neutron_count(10, 5), 0);
    }

    #[test]
    fn binding_deficit_helium4() {
        // He-4: mass ≈ 4.0026 u, Z=2, A=4 → ~28 MeV binding
        let be = mass_to_binding_deficit_mev(2, 4, 4.002_602);
        assert!(be > 25.0 && be < 35.0, "He-4 binding: {be} MeV");
    }

    #[test]
    fn decay_constant_basic() {
        // t½ = 1 s → λ = ln(2) ≈ 0.693
        let lam = half_life_to_decay_constant(1.0);
        assert!((lam - core::f64::consts::LN_2).abs() < 1e-10);
    }

    #[test]
    fn decay_constant_zero() {
        assert!(half_life_to_decay_constant(0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn transitions_to_wavelengths_basic() {
        let energies = [10.2, 12.09, 1.89]; // Lyman-α, Lyman-β, H-α
        let wl = transitions_to_wavelengths(&energies);
        assert_eq!(wl.len(), 3);
        assert!((wl[0] - 121.6).abs() < 1.0);
        assert!((wl[2] - 656.0).abs() < 5.0);
    }

    #[test]
    fn transitions_skip_negative() {
        let wl = transitions_to_wavelengths(&[-1.0, 0.0, 10.2]);
        assert_eq!(wl.len(), 1);
    }

    #[test]
    fn hyperfine_hydrogen() {
        // Hydrogen: I = 1/2 → 2 states
        let s = nuclear_spin_to_hyperfine_scale(0.5);
        assert!((s - 2.0).abs() < 0.01);
    }

    #[test]
    fn hyperfine_zero_spin() {
        assert!(nuclear_spin_to_hyperfine_scale(0.0).abs() < f64::EPSILON);
    }
}
