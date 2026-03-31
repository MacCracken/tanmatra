//! Cross-crate bridges — convert primitive values from other AGNOS science crates
//! into tanmatra atomic/nuclear physics parameters and vice versa.
//!
//! Always available — takes primitive values (f64), no science crate deps.

use alloc::vec::Vec;
use serde::{Deserialize, Serialize};

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

// ── Jyotish bridges (astronomy/time) ──────────────────────────────────────

/// Convert TAI seconds to TT (Terrestrial Time) seconds.
///
/// TT = TAI + 32.184 s (IAU 1991). This bridge allows jyotish to replace
/// polynomial ΔT approximations with exact TAI↔TT conversion.
#[must_use]
#[inline]
pub fn tai_to_tt(tai_seconds: f64) -> f64 {
    crate::timekeeping::tai_to_tt(tai_seconds)
}

// ── Chrono bridge (datetime interop) ─────────────────────────────────────

/// Convert a UTC date to the TAI−UTC offset (leap seconds) as integer seconds.
///
/// Allows chrono consumers to compute the TAI instant corresponding to a
/// `chrono::DateTime<Utc>` by adding this offset to the UTC epoch seconds.
///
/// `year`: Gregorian year. `month`: 1–12.
/// Returns ΔAT = TAI − UTC in whole seconds (0 before 1972).
#[must_use]
#[inline]
pub fn utc_date_to_tai_offset(year: i32, month: u8) -> i32 {
    crate::timekeeping::leap_seconds_at(year, month)
}

/// Convert TAI seconds-since-epoch to UTC seconds-since-epoch.
///
/// UTC = TAI − ΔAT. The caller supplies the UTC year/month for leap-second
/// lookup (chrono can derive these from the rough UTC estimate).
#[must_use]
#[inline]
pub fn tai_to_utc_seconds(tai_seconds: f64, utc_year: i32, utc_month: u8) -> f64 {
    tai_seconds - crate::timekeeping::leap_seconds_at(utc_year, utc_month) as f64
}

/// Convert UTC seconds-since-epoch to TAI seconds-since-epoch.
///
/// TAI = UTC + ΔAT.
#[must_use]
#[inline]
pub fn utc_to_tai_seconds(utc_seconds: f64, year: i32, month: u8) -> f64 {
    utc_seconds + crate::timekeeping::leap_seconds_at(year, month) as f64
}

// ── Hisab-mimamsa bridge (relativity/gravitation) ────────────────────────

/// Gravitational time dilation factor for a clock at `orbital_radius_m`
/// relative to a ground clock at Earth's surface.
///
/// Returns the fractional rate difference Δf/f (positive = satellite clock
/// runs faster than ground). Wraps the Schwarzschild correction from
/// [`crate::timekeeping`].
///
/// For GPS orbit: returns ≈ +4.465e-10 (≈ +38.6 μs/day).
#[must_use]
#[inline]
pub fn gravitational_time_dilation(orbital_radius_m: f64, orbital_velocity_m_s: f64) -> f64 {
    // Convert μs/day → dimensionless fractional rate
    let us_per_day = crate::timekeeping::schwarzschild_clock_correction_us_per_day(
        orbital_radius_m,
        orbital_velocity_m_s,
    );
    us_per_day / (86_400.0 * 1e6)
}

/// Gravitational time dilation accumulated over a duration (seconds).
///
/// Returns the clock offset in seconds: Δt = (Δf/f) × duration_s.
/// Positive means the orbiting clock gains time relative to ground.
#[must_use]
#[inline]
pub fn gravitational_time_offset_s(
    orbital_radius_m: f64,
    orbital_velocity_m_s: f64,
    duration_s: f64,
) -> f64 {
    gravitational_time_dilation(orbital_radius_m, orbital_velocity_m_s) * duration_s
}

// ── Falak bridge (orbital mechanics/ephemeris) ───────────────────────────

/// Convert TAI seconds since the TAI epoch (1958-01-01) to Julian Date (TT).
///
/// JD(TT) = TAI_EPOCH_JD_TT + tai_seconds / 86400, where the epoch constant
/// already encodes the 32.184 s TT−TAI offset.
#[must_use]
#[inline]
pub fn tai_seconds_to_jd_tt(tai_seconds: f64) -> f64 {
    // 1958-01-01T00:00:00 TAI in JD(TT) = 2436204.5 + 32.184/86400
    const TAI_EPOCH_JD_TT: f64 = 2_436_204.500_372_5;
    TAI_EPOCH_JD_TT + tai_seconds / 86_400.0
}

/// Convert Julian Date (TT) to TAI seconds since the TAI epoch (1958-01-01).
///
/// Inverse of [`tai_seconds_to_jd_tt`].
#[must_use]
#[inline]
pub fn jd_tt_to_tai_seconds(jd_tt: f64) -> f64 {
    const TAI_EPOCH_JD_TT: f64 = 2_436_204.500_372_5;
    (jd_tt - TAI_EPOCH_JD_TT) * 86_400.0
}

/// Convert TAI seconds to Modified Julian Date (TT).
///
/// MJD = JD − 2400000.5.
#[must_use]
#[inline]
pub fn tai_seconds_to_mjd_tt(tai_seconds: f64) -> f64 {
    tai_seconds_to_jd_tt(tai_seconds) - 2_400_000.5
}

// ── Kiran/Joshua bridge (simulation time) ────────────────────────────────

/// A simulation clock anchored to physical atomic time.
///
/// Wraps a TAI-referenced origin with a time-scale multiplier and
/// accumulated offset, supporting pause, fast-forward, and slow-motion.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SimulationClock {
    /// TAI seconds at which the clock was started/last-resumed.
    origin_tai_s: f64,
    /// Accumulated simulation seconds at the point the clock was last set.
    accumulated_s: f64,
    /// Simulation time multiplier (1.0 = real-time, 0.0 = paused, 2.0 = 2× speed).
    multiplier: f64,
}

impl SimulationClock {
    /// Creates a new simulation clock anchored at `origin_tai_s` with 1× speed.
    #[must_use]
    pub fn new(origin_tai_s: f64) -> Self {
        Self {
            origin_tai_s,
            accumulated_s: 0.0,
            multiplier: 1.0,
        }
    }

    /// Returns the simulation time (seconds) at a given TAI instant.
    ///
    /// sim_time = accumulated + multiplier × (current_tai − origin_tai).
    #[must_use]
    #[inline]
    pub fn simulation_time(&self, current_tai_s: f64) -> f64 {
        self.accumulated_s + self.multiplier * (current_tai_s - self.origin_tai_s)
    }

    /// Returns the current multiplier.
    #[must_use]
    #[inline]
    pub fn multiplier(&self) -> f64 {
        self.multiplier
    }

    /// Returns `true` if the clock is paused (multiplier == 0).
    #[must_use]
    #[inline]
    pub fn is_paused(&self) -> bool {
        self.multiplier == 0.0
    }

    /// Sets the time-scale multiplier, freezing the accumulated time at the
    /// given TAI instant before changing the rate.
    ///
    /// Use `0.0` to pause, `1.0` for real-time, `>1.0` for fast-forward.
    #[must_use]
    pub fn set_multiplier(self, current_tai_s: f64, new_multiplier: f64) -> Self {
        let acc = self.simulation_time(current_tai_s);
        Self {
            origin_tai_s: current_tai_s,
            accumulated_s: acc,
            multiplier: new_multiplier,
        }
    }

    /// Pauses the clock at the given TAI instant.
    #[must_use]
    #[inline]
    pub fn pause(self, current_tai_s: f64) -> Self {
        self.set_multiplier(current_tai_s, 0.0)
    }

    /// Resumes the clock at 1× speed from the given TAI instant.
    #[must_use]
    #[inline]
    pub fn resume(self, current_tai_s: f64) -> Self {
        self.set_multiplier(current_tai_s, 1.0)
    }
}

// ── Bhava bridge (circadian/rhythm time context) ─────────────────────────

/// Time context for circadian, rhythm, and growth modules.
///
/// Distinguishes between real-time (wall-clock), simulated (game time
/// with a potentially different rate), and paused states so that
/// biological/rhythmic models can adapt their behaviour accordingly.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum TimeContext {
    /// Wall-clock time: simulation runs at 1× real-time, anchored to TAI.
    ///
    /// `tai_seconds`: current TAI seconds since epoch.
    RealTime {
        /// Current TAI seconds since 1958-01-01T00:00:00 TAI.
        tai_seconds: f64,
    },
    /// Simulated time: simulation runs at an arbitrary rate.
    ///
    /// Circadian modules should use `sim_elapsed_s` for phase accumulation
    /// rather than wall-clock time.
    Simulated {
        /// Current TAI seconds (wall-clock reference).
        tai_seconds: f64,
        /// Elapsed simulation seconds since simulation start.
        sim_elapsed_s: f64,
        /// Simulation speed multiplier (e.g., 2.0 = double speed).
        multiplier: f64,
    },
    /// Paused: no time passes for rhythmic/growth processes.
    ///
    /// Modules receiving `Paused` should freeze their state.
    Paused {
        /// TAI seconds at which the pause began.
        tai_seconds: f64,
        /// Frozen simulation seconds.
        sim_elapsed_s: f64,
    },
}

impl TimeContext {
    /// Creates a real-time context from a TAI instant.
    #[must_use]
    #[inline]
    pub fn real_time(tai_seconds: f64) -> Self {
        Self::RealTime { tai_seconds }
    }

    /// Creates a simulated context from a [`SimulationClock`] and current TAI.
    #[must_use]
    pub fn from_simulation_clock(clock: &SimulationClock, current_tai_s: f64) -> Self {
        if clock.is_paused() {
            Self::Paused {
                tai_seconds: current_tai_s,
                sim_elapsed_s: clock.simulation_time(current_tai_s),
            }
        } else {
            Self::Simulated {
                tai_seconds: current_tai_s,
                sim_elapsed_s: clock.simulation_time(current_tai_s),
                multiplier: clock.multiplier(),
            }
        }
    }

    /// Returns the TAI seconds regardless of context variant.
    #[must_use]
    pub fn tai_seconds(&self) -> f64 {
        match *self {
            Self::RealTime { tai_seconds }
            | Self::Simulated { tai_seconds, .. }
            | Self::Paused { tai_seconds, .. } => tai_seconds,
        }
    }

    /// Returns the effective simulation elapsed time.
    ///
    /// For `RealTime`, this equals the TAI seconds (simulation ≡ real).
    /// For `Simulated` and `Paused`, returns the simulation elapsed time.
    #[must_use]
    pub fn effective_elapsed_s(&self) -> f64 {
        match *self {
            Self::RealTime { tai_seconds } => tai_seconds,
            Self::Simulated { sim_elapsed_s, .. } | Self::Paused { sim_elapsed_s, .. } => {
                sim_elapsed_s
            }
        }
    }

    /// Returns `true` if the context is paused.
    #[must_use]
    pub fn is_paused(&self) -> bool {
        matches!(self, Self::Paused { .. })
    }

    /// Returns `true` if running in real-time (1× speed).
    #[must_use]
    pub fn is_real_time(&self) -> bool {
        matches!(self, Self::RealTime { .. })
    }

    /// Returns the simulation multiplier (1.0 for real-time, 0.0 for paused).
    #[must_use]
    pub fn multiplier(&self) -> f64 {
        match *self {
            Self::RealTime { .. } => 1.0,
            Self::Simulated { multiplier, .. } => multiplier,
            Self::Paused { .. } => 0.0,
        }
    }
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

    // ── Chrono bridge tests ──────────────────────────────────────────────

    #[test]
    fn utc_date_to_tai_offset_2020() {
        // 2020: ΔAT = 37 s (last leap second was 2017-01-01)
        assert_eq!(utc_date_to_tai_offset(2020, 6), 37);
    }

    #[test]
    fn utc_date_to_tai_offset_before_1972() {
        assert_eq!(utc_date_to_tai_offset(1960, 1), 0);
    }

    #[test]
    fn tai_utc_roundtrip() {
        let tai = 1_000_000.0;
        let utc = tai_to_utc_seconds(tai, 2020, 1);
        let back = utc_to_tai_seconds(utc, 2020, 1);
        assert!((back - tai).abs() < 1e-10);
    }

    #[test]
    fn utc_tai_offset_matches_timekeeping() {
        let offset = utc_date_to_tai_offset(2015, 8);
        assert_eq!(offset, 36); // After 2015-07-01 leap second
    }

    // ── Hisab-mimamsa bridge tests ──────────────────────────────────────

    #[test]
    fn gravitational_dilation_gps() {
        // GPS orbit: R ≈ 26,560 km, v ≈ 3874 m/s → ~4.5e-10 fractional
        let frac = gravitational_time_dilation(26_560_000.0, 3874.0);
        assert!(frac > 4e-10 && frac < 5e-10, "frac = {frac}");
    }

    #[test]
    fn gravitational_offset_one_day() {
        let offset = gravitational_time_offset_s(26_560_000.0, 3874.0, 86_400.0);
        // ~38.6 μs/day = ~3.86e-5 s
        assert!((offset - 3.86e-5).abs() < 1e-5, "offset = {offset} s");
    }

    #[test]
    fn gravitational_dilation_ground() {
        // At surface radius with zero velocity → should be ~0 (reference)
        let frac = gravitational_time_dilation(6_371_000.0, 0.0);
        assert!(frac.abs() < 1e-12, "ground frac = {frac}");
    }

    // ── Falak bridge tests ──────────────────────────────────────────────

    #[test]
    fn tai_to_jd_roundtrip() {
        let tai = 500_000_000.0;
        let jd = tai_seconds_to_jd_tt(tai);
        let back = jd_tt_to_tai_seconds(jd);
        // f64 JD precision: ~2.4M days × 86400 s/day ≈ 2e11, so ~1e-4 s roundtrip error
        assert!((back - tai).abs() < 1e-3, "roundtrip error: {}", back - tai);
    }

    #[test]
    fn tai_epoch_jd_tt() {
        // At TAI = 0, JD(TT) should be ~2436204.5 + 32.184/86400
        let jd = tai_seconds_to_jd_tt(0.0);
        assert!((jd - 2_436_204.500_372_5).abs() < 1e-6, "epoch JD = {jd}");
    }

    #[test]
    fn mjd_offset() {
        let tai = 0.0;
        let jd = tai_seconds_to_jd_tt(tai);
        let mjd = tai_seconds_to_mjd_tt(tai);
        assert!((mjd - (jd - 2_400_000.5)).abs() < 1e-10);
    }

    #[test]
    fn j2000_epoch_check() {
        // J2000.0 = JD 2451545.0 (TT). TAI seconds for that:
        let j2000_tai = jd_tt_to_tai_seconds(2_451_545.0);
        // Should be ~1325376000 - 32.184 ≈ positive large number (post-1958)
        assert!(j2000_tai > 0.0, "J2000 TAI = {j2000_tai}");
        // Roundtrip
        let jd_back = tai_seconds_to_jd_tt(j2000_tai);
        assert!((jd_back - 2_451_545.0).abs() < 1e-6);
    }

    // ── Kiran/Joshua bridge tests ───────────────────────────────────────

    #[test]
    fn simulation_clock_realtime() {
        let clock = SimulationClock::new(100.0);
        assert!((clock.simulation_time(200.0) - 100.0).abs() < 1e-10);
        assert!((clock.multiplier() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn simulation_clock_fast_forward() {
        let clock = SimulationClock::new(0.0).set_multiplier(0.0, 10.0);
        // At TAI = 5.0, sim time = 0 + 10 × 5 = 50
        assert!((clock.simulation_time(5.0) - 50.0).abs() < 1e-10);
    }

    #[test]
    fn simulation_clock_pause_resume() {
        let clock = SimulationClock::new(0.0);
        // Run for 10s at 1×, then pause at TAI=10
        let paused = clock.pause(10.0);
        assert!(paused.is_paused());
        // Sim time frozen at 10.0 regardless of wall-clock advance
        assert!((paused.simulation_time(100.0) - 10.0).abs() < 1e-10);
        // Resume at TAI=20
        let resumed = paused.resume(20.0);
        assert!(!resumed.is_paused());
        // At TAI=25, 5s of real-time passed since resume → sim = 10 + 5 = 15
        assert!((resumed.simulation_time(25.0) - 15.0).abs() < 1e-10);
    }

    #[test]
    fn simulation_clock_slow_motion() {
        let clock = SimulationClock::new(0.0).set_multiplier(0.0, 0.5);
        // At TAI = 10, sim = 0 + 0.5 × 10 = 5
        assert!((clock.simulation_time(10.0) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn serde_roundtrip_simulation_clock() {
        let clock = SimulationClock::new(1000.0).set_multiplier(1000.0, 2.0);
        let json = serde_json::to_string(&clock).unwrap();
        let back: SimulationClock = serde_json::from_str(&json).unwrap();
        assert_eq!(clock, back);
    }

    // ── Bhava bridge tests ──────────────────────────────────────────────

    #[test]
    fn time_context_real_time() {
        let ctx = TimeContext::real_time(1_000_000.0);
        assert!(ctx.is_real_time());
        assert!(!ctx.is_paused());
        assert!((ctx.multiplier() - 1.0).abs() < 1e-10);
        assert!((ctx.tai_seconds() - 1_000_000.0).abs() < 1e-10);
    }

    #[test]
    fn time_context_from_clock_running() {
        let clock = SimulationClock::new(0.0).set_multiplier(0.0, 3.0);
        let ctx = TimeContext::from_simulation_clock(&clock, 10.0);
        assert!(!ctx.is_paused());
        assert!(!ctx.is_real_time());
        assert!((ctx.multiplier() - 3.0).abs() < 1e-10);
        assert!((ctx.effective_elapsed_s() - 30.0).abs() < 1e-10);
    }

    #[test]
    fn time_context_from_clock_paused() {
        let clock = SimulationClock::new(0.0).pause(5.0);
        let ctx = TimeContext::from_simulation_clock(&clock, 100.0);
        assert!(ctx.is_paused());
        assert!((ctx.multiplier()).abs() < 1e-10);
        assert!((ctx.effective_elapsed_s() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn time_context_tai_all_variants() {
        let tai = 42.0;
        let rt = TimeContext::RealTime { tai_seconds: tai };
        let sim = TimeContext::Simulated {
            tai_seconds: tai,
            sim_elapsed_s: 100.0,
            multiplier: 2.0,
        };
        let paused = TimeContext::Paused {
            tai_seconds: tai,
            sim_elapsed_s: 50.0,
        };
        assert!((rt.tai_seconds() - tai).abs() < 1e-10);
        assert!((sim.tai_seconds() - tai).abs() < 1e-10);
        assert!((paused.tai_seconds() - tai).abs() < 1e-10);
    }

    #[test]
    fn serde_roundtrip_time_context() {
        let ctx = TimeContext::Simulated {
            tai_seconds: 1000.0,
            sim_elapsed_s: 5000.0,
            multiplier: 2.5,
        };
        let json = serde_json::to_string(&ctx).unwrap();
        let back: TimeContext = serde_json::from_str(&json).unwrap();
        assert_eq!(ctx, back);
    }

    #[test]
    fn serde_roundtrip_time_context_paused() {
        let ctx = TimeContext::Paused {
            tai_seconds: 500.0,
            sim_elapsed_s: 200.0,
        };
        let json = serde_json::to_string(&ctx).unwrap();
        let back: TimeContext = serde_json::from_str(&json).unwrap();
        assert_eq!(ctx, back);
    }
}
