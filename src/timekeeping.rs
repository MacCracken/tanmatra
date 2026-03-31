//! Frequency standards, atomic time scales, and relativistic clock corrections.
//!
//! This module provides:
//! - [`FrequencyStandard`](crate::timekeeping::FrequencyStandard): Atomic frequency standards (Cs, Rb, H, Sr, Yb)
//! - [`TimeScale`](crate::timekeeping::TimeScale): Atomic and astronomical time scales (TAI, UTC, TT, GPS, TCB, TCG)
//! - [`AtomicInstant`](crate::timekeeping::AtomicInstant): TAI-referenced instant with sub-nanosecond precision
//! - Leap second table (IERS Bulletin C, 1972–2017)
//! - Relativistic clock corrections (gravitational redshift, Sagnac, Doppler)
//!
//! ## Data Sources
//!
//! - **BIPM/CIPM 2021**: Secondary representations of the SI second
//! - **CODATA 2022**: Fundamental constants
//! - **IERS Bulletin C**: Leap second announcements
//! - **IAU 2000/2006**: Time scale definitions (L_G, L_C)

use crate::constants::{C, EARTH_ROTATION_RAD_S, GM_EARTH, STANDARD_GRAVITY};
use serde::{Deserialize, Serialize};

// ── Frequency Standards ─────────────────────────────────────────────────────

/// Atomic frequency standards used to define or approximate the SI second.
///
/// Each variant represents a real atomic transition used in precision
/// timekeeping. The cesium-133 hyperfine transition defines the SI second
/// (CGPM 1967, reaffirmed 2019). The optical standards are BIPM 2021
/// secondary representations of the second.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum FrequencyStandard {
    /// Cesium-133 hyperfine transition (SI second definition, CGPM 1967/2019).
    ///
    /// Transition: ⁶S₁/₂ (F=3 → F=4), frequency = 9 192 631 770 Hz (exact).
    Cesium133,
    /// Rubidium-87 hyperfine transition.
    ///
    /// Transition: ⁵S₁/₂ (F=1 → F=2), frequency ≈ 6 834 682 610.904 312 Hz.
    /// Reference: BIPM recommended value.
    Rubidium87,
    /// Hydrogen maser (1420 MHz, 21 cm line).
    ///
    /// Transition: 1S₁/₂ (F=0 → F=1), frequency ≈ 1 420 405 751.768 Hz.
    /// Reference: NIST.
    HydrogenMaser,
    /// Strontium-87 optical lattice clock (BIPM 2021 secondary representation).
    ///
    /// Transition: ¹S₀ → ³P₀, frequency ≈ 429 228 004 229 873.2 Hz.
    /// Reference: BIPM CIPM 2021.
    StrontiumOptical,
    /// Ytterbium-171 optical lattice clock (BIPM 2021 secondary representation).
    ///
    /// Transition: ¹S₀ → ³P₀, frequency ≈ 518 295 836 590 863.6 Hz.
    /// Reference: BIPM CIPM 2021.
    YtterbiumOptical,
}

impl FrequencyStandard {
    /// Transition frequency in Hz.
    ///
    /// For Cs-133, this is exact by definition (SI second).
    /// Other values are from BIPM 2021 recommended values.
    #[must_use]
    #[inline]
    pub fn transition_frequency_hz(self) -> f64 {
        match self {
            Self::Cesium133 => 9_192_631_770.0,
            Self::Rubidium87 => 6_834_682_610.904_312,
            Self::HydrogenMaser => 1_420_405_751.768,
            Self::StrontiumOptical => 429_228_004_229_873.2,
            Self::YtterbiumOptical => 518_295_836_590_863.6,
        }
    }

    /// Transition wavelength in meters (λ = c / f).
    #[must_use]
    #[inline]
    pub fn transition_wavelength_m(self) -> f64 {
        C / self.transition_frequency_hz()
    }

    /// Approximate quality factor Q = f × τ_coherence.
    ///
    /// These are order-of-magnitude values representative of each technology.
    #[must_use]
    #[inline]
    pub fn quality_factor(self) -> f64 {
        match self {
            Self::Cesium133 => 1e10,
            Self::Rubidium87 => 1e10,
            Self::HydrogenMaser => 1e9,
            Self::StrontiumOptical => 1e17,
            Self::YtterbiumOptical => 1e17,
        }
    }

    /// Typical fractional stability (Allan deviation at 1 second).
    ///
    /// These are representative values for each clock technology.
    #[must_use]
    #[inline]
    pub fn fractional_stability(self) -> f64 {
        match self {
            Self::Cesium133 => 1e-13,
            Self::Rubidium87 => 1e-13,
            Self::HydrogenMaser => 1e-15,
            Self::StrontiumOptical => 1e-18,
            Self::YtterbiumOptical => 1e-18,
        }
    }
}

/// Returns a static slice of all frequency standards.
#[must_use]
pub fn all_frequency_standards() -> &'static [FrequencyStandard] {
    &[
        FrequencyStandard::Cesium133,
        FrequencyStandard::Rubidium87,
        FrequencyStandard::HydrogenMaser,
        FrequencyStandard::StrontiumOptical,
        FrequencyStandard::YtterbiumOptical,
    ]
}

// ── Atomic Time Scales ──────────────────────────────────────────────────────

/// Atomic and astronomical time scales.
///
/// These represent the standard time scales used in precision timekeeping,
/// geodesy, and astronomy. Conversions between them are defined by
/// IAU resolutions and BIPM conventions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum TimeScale {
    /// International Atomic Time (Temps Atomique International).
    TAI,
    /// Coordinated Universal Time (with leap seconds).
    UTC,
    /// Terrestrial Time (successor to Ephemeris Time, IAU 1991).
    TT,
    /// GPS Time (epoch: 1980-01-06T00:00:00 UTC).
    GPS,
    /// Barycentric Coordinate Time (IAU 2006).
    TCB,
    /// Geocentric Coordinate Time (IAU 2000).
    TCG,
}

/// TT − TAI offset in seconds (IAU 1991, exact).
pub const TAI_TT_OFFSET_S: f64 = 32.184;

/// TAI − GPS offset in seconds (exact).
pub const TAI_GPS_OFFSET_S: f64 = 19.0;

/// TCG/TT rate difference L_G (IAU 2000 Resolution B1.9, defining constant).
///
/// dTCG/dTT = 1 + L_G.
pub const LG_RATE: f64 = 6.969_290_134e-10;

/// TCB/TCG rate difference L_C (IAU 2006 Resolution B3).
///
/// dTCB/dTCG ≈ 1 + L_C.
pub const LC_RATE: f64 = 1.480_826_867_41e-8;

/// Convert TAI seconds to TT seconds.
///
/// TT = TAI + 32.184 s (IAU 1991).
#[must_use]
#[inline]
pub fn tai_to_tt(tai_seconds: f64) -> f64 {
    tai_seconds + TAI_TT_OFFSET_S
}

/// Convert TT seconds to TAI seconds.
///
/// TAI = TT − 32.184 s.
#[must_use]
#[inline]
pub fn tt_to_tai(tt_seconds: f64) -> f64 {
    tt_seconds - TAI_TT_OFFSET_S
}

/// Convert TAI seconds to GPS seconds.
///
/// GPS = TAI − 19 s.
#[must_use]
#[inline]
pub fn tai_to_gps(tai_seconds: f64) -> f64 {
    tai_seconds - TAI_GPS_OFFSET_S
}

/// Convert GPS seconds to TAI seconds.
///
/// TAI = GPS + 19 s.
#[must_use]
#[inline]
pub fn gps_to_tai(gps_seconds: f64) -> f64 {
    gps_seconds + TAI_GPS_OFFSET_S
}

// ── Leap Second Table ───────────────────────────────────────────────────────

/// Leap second entries from IERS Bulletin C.
///
/// Format: (year, month, delta_at_after)
/// - month: 1 = January, 7 = July
/// - delta_at_after: TAI − UTC after the leap second (in whole seconds)
///
/// Source: IERS Bulletin C, complete through 2017-01-01.
/// As of IERS Bulletin C 69 (2025), no further leap seconds have been
/// announced; TAI − UTC remains 37 s.
const LEAP_SECONDS: &[(i32, u8, i32)] = &[
    (1972, 1, 10),
    (1972, 7, 11),
    (1973, 1, 12),
    (1974, 1, 13),
    (1975, 1, 14),
    (1976, 1, 15),
    (1977, 1, 16),
    (1978, 1, 17),
    (1979, 1, 18),
    (1980, 1, 19),
    (1981, 7, 20),
    (1982, 7, 21),
    (1983, 7, 22),
    (1985, 7, 23),
    (1988, 1, 24),
    (1990, 1, 25),
    (1991, 1, 26),
    (1992, 7, 27),
    (1993, 7, 28),
    (1994, 7, 29),
    (1996, 1, 30),
    (1997, 7, 31),
    (1999, 1, 32),
    (2006, 1, 33),
    (2009, 1, 34),
    (2012, 7, 35),
    (2015, 7, 36),
    (2017, 1, 37),
];

/// Returns TAI − UTC (ΔAT) in whole seconds for a given date.
///
/// Scans the leap second table to find the most recent entry at or before
/// the given (year, month). Returns 0 for dates before 1972-01-01.
///
/// # Arguments
///
/// * `year` — Calendar year (e.g. 2024).
/// * `month` — Month number (1–12).
#[must_use]
pub fn leap_seconds_at(year: i32, month: u8) -> i32 {
    let mut result = 0;
    for &(y, m, delta_at) in LEAP_SECONDS {
        if year > y || (year == y && month >= m) {
            result = delta_at;
        } else {
            break;
        }
    }
    result
}

/// Returns TAI − UTC offset as f64 for a given date.
///
/// Equivalent to [`leap_seconds_at`] but as floating-point seconds.
#[must_use]
#[inline]
pub fn tai_to_utc_offset(year: i32, month: u8) -> f64 {
    leap_seconds_at(year, month) as f64
}

/// Returns UTC − TAI offset as f64 for a given date (negated).
///
/// This is the additive correction to convert TAI to UTC:
/// UTC = TAI − ΔAT, so the offset is −ΔAT.
#[must_use]
#[inline]
pub fn utc_to_tai_offset(year: i32, month: u8) -> f64 {
    -(leap_seconds_at(year, month) as f64)
}

// ── AtomicInstant ───────────────────────────────────────────────────────────

/// A TAI-referenced instant with sub-nanosecond precision.
///
/// Internally stores whole seconds (i64) and fractional nanoseconds (u32).
/// Epoch: 1958-01-01T00:00:00 TAI (the conventional TAI epoch).
///
/// # Examples
///
/// ```
/// # // Note: this doctest won't run since module isn't pub yet
/// use tanmatra::timekeeping::AtomicInstant;
/// let t = AtomicInstant::new(1_000_000_000, 500_000_000);
/// assert!(t.tai_seconds() > 999_999_999.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct AtomicInstant {
    /// Whole seconds since 1958-01-01T00:00:00 TAI.
    seconds: i64,
    /// Sub-second nanoseconds \[0, 999_999_999\].
    nanos: u32,
}

/// Maximum valid nanosecond value.
const MAX_NANOS: u32 = 999_999_999;

/// Nanoseconds per second.
const NANOS_PER_SEC: f64 = 1_000_000_000.0;

impl AtomicInstant {
    /// Creates a new `AtomicInstant`.
    ///
    /// Clamps `nanos` to \[0, 999_999_999\]. If nanos exceeds this range,
    /// the excess is carried into whole seconds.
    #[must_use]
    pub fn new(seconds: i64, nanos: u32) -> Self {
        if nanos > MAX_NANOS {
            let extra_secs = (nanos / (MAX_NANOS + 1)) as i64;
            let remaining = nanos % (MAX_NANOS + 1);
            Self {
                seconds: seconds.saturating_add(extra_secs),
                nanos: remaining,
            }
        } else {
            Self { seconds, nanos }
        }
    }

    /// Creates an `AtomicInstant` from floating-point TAI seconds since epoch.
    #[must_use]
    pub fn from_tai_seconds(s: f64) -> Self {
        let whole = libm::floor(s) as i64;
        let frac = s - libm::floor(s);
        // frac is in [0.0, 1.0) after floor subtraction, so the cast is safe.
        let nanos_f = frac * NANOS_PER_SEC;
        #[allow(clippy::cast_sign_loss)]
        let nanos = if nanos_f < 0.0 { 0u32 } else { nanos_f as u32 };
        Self::new(whole, nanos)
    }

    /// Returns the TAI time as floating-point seconds since epoch.
    #[must_use]
    pub fn tai_seconds(&self) -> f64 {
        self.seconds as f64 + self.nanos as f64 / NANOS_PER_SEC
    }

    /// Returns the time as TT seconds since epoch.
    ///
    /// TT = TAI + 32.184 s.
    #[must_use]
    #[inline]
    pub fn to_tt_seconds(&self) -> f64 {
        tai_to_tt(self.tai_seconds())
    }

    /// Returns the time as GPS seconds since epoch.
    ///
    /// GPS = TAI − 19 s.
    #[must_use]
    #[inline]
    pub fn to_gps_seconds(&self) -> f64 {
        tai_to_gps(self.tai_seconds())
    }

    /// Returns the duration in seconds between this instant and an earlier one.
    ///
    /// If `earlier` is actually later than `self`, the result will be negative.
    #[must_use]
    pub fn duration_since(&self, earlier: &Self) -> f64 {
        let dsec = self.seconds - earlier.seconds;
        let dnanos = self.nanos as i64 - earlier.nanos as i64;
        dsec as f64 + dnanos as f64 / NANOS_PER_SEC
    }

    /// Returns a new instant advanced by `dt` seconds.
    #[must_use]
    pub fn add_seconds(&self, dt: f64) -> Self {
        Self::from_tai_seconds(self.tai_seconds() + dt)
    }
}

// ── Relativistic Clock Corrections ──────────────────────────────────────────

/// Gravitational redshift fractional frequency shift.
///
/// Δf/f = −g Δh / c², where g = 9.80665 m/s² (standard gravity).
///
/// A positive `delta_h_m` (receiver higher than emitter) gives a *negative*
/// fractional shift (blueshift when received from below).
///
/// # Arguments
///
/// * `delta_h_m` — Height difference in meters (positive = upward).
#[must_use]
#[inline]
pub fn gravitational_redshift(delta_h_m: f64) -> f64 {
    -STANDARD_GRAVITY * delta_h_m / (C * C)
}

/// Mean Earth radius in meters (IERS 2010).
const EARTH_MEAN_RADIUS_M: f64 = 6.371_000e6;

/// Total relativistic clock correction for an orbiting satellite relative
/// to a ground clock.
///
/// Combines gravitational blueshift (satellite clock runs faster at higher
/// gravitational potential) and velocity time dilation (moving clock runs
/// slower):
///
/// - Gravitational: +GM/c² × (1/R_earth − 1/R_orbit) × 86400 × 10⁶ μs/day
/// - Velocity: −v²/(2c²) × 86400 × 10⁶ μs/day
///
/// Uses GM_earth = 3.986004418 × 10¹⁴ m³/s² (IERS 2010).
///
/// # Arguments
///
/// * `orbital_radius_m` — Orbital radius in meters (from Earth center).
/// * `orbital_velocity_m_s` — Orbital velocity in m/s.
///
/// # Returns
///
/// Net clock correction in microseconds per day. Positive means the
/// satellite clock runs *fast* relative to a ground clock.
///
/// # Example: GPS satellite
///
/// R ≈ 26 560 km, v ≈ 3874 m/s → gravitational +45.85 μs/day,
/// velocity −7.21 μs/day, net ≈ +38.6 μs/day.
#[must_use]
pub fn schwarzschild_clock_correction_us_per_day(
    orbital_radius_m: f64,
    orbital_velocity_m_s: f64,
) -> f64 {
    let c2 = C * C;
    let seconds_per_day: f64 = 86400.0;
    let us_per_s: f64 = 1e6;

    // Gravitational potential difference between ground and orbit
    let grav = GM_EARTH / c2 * (1.0 / EARTH_MEAN_RADIUS_M - 1.0 / orbital_radius_m);
    let vel = orbital_velocity_m_s * orbital_velocity_m_s / (2.0 * c2);

    (grav - vel) * seconds_per_day * us_per_s
}

/// Second-order (transverse) Doppler shift.
///
/// Δf/f = −β²/2, where β = v/c.
///
/// This is the purely relativistic time dilation effect for a
/// moving clock (no radial component needed).
///
/// # Arguments
///
/// * `beta` — v/c (dimensionless speed parameter).
#[must_use]
#[inline]
pub fn second_order_doppler_shift(beta: f64) -> f64 {
    -beta * beta / 2.0
}

/// Sagnac effect correction for a signal traversing a path on a rotating Earth.
///
/// Δt = 4 Ω A sin(φ) / c², where:
/// - Ω = 7.2921150 × 10⁻⁵ rad/s (Earth rotation rate)
/// - A = projected area enclosed by the signal path (m²)
/// - φ = geodetic latitude (rad)
///
/// # Arguments
///
/// * `latitude_rad` — Geodetic latitude in radians.
/// * `area_m2` — Area enclosed by the signal path in m².
///
/// # Returns
///
/// Sagnac correction in nanoseconds.
#[must_use]
pub fn sagnac_correction_ns(latitude_rad: f64, area_m2: f64) -> f64 {
    let c2 = C * C;
    let dt_s = 4.0 * EARTH_ROTATION_RAD_S * area_m2 * libm::sin(latitude_rad) / c2;
    dt_s * 1e9
}

// ── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── FrequencyStandard ───────────────────────────────────────────────

    #[test]
    fn cesium_frequency_exact() {
        assert!(
            (FrequencyStandard::Cesium133.transition_frequency_hz() - 9_192_631_770.0).abs()
                < f64::EPSILON
        );
    }

    #[test]
    fn rubidium_frequency() {
        let f = FrequencyStandard::Rubidium87.transition_frequency_hz();
        assert!((f - 6_834_682_610.904_312).abs() < 1.0);
    }

    #[test]
    fn hydrogen_maser_frequency() {
        let f = FrequencyStandard::HydrogenMaser.transition_frequency_hz();
        assert!((f - 1_420_405_751.768).abs() < 1.0);
    }

    #[test]
    fn strontium_optical_frequency() {
        let f = FrequencyStandard::StrontiumOptical.transition_frequency_hz();
        assert!((f - 429_228_004_229_873.2).abs() < 10.0);
    }

    #[test]
    fn ytterbium_optical_frequency() {
        let f = FrequencyStandard::YtterbiumOptical.transition_frequency_hz();
        assert!((f - 518_295_836_590_863.6).abs() < 10.0);
    }

    #[test]
    fn cesium_wavelength() {
        // Cs-133: λ = c / 9.19e9 ≈ 3.26 cm (microwave)
        let wl = FrequencyStandard::Cesium133.transition_wavelength_m();
        assert!((wl - 0.0326).abs() < 0.001);
    }

    #[test]
    fn strontium_wavelength() {
        // Sr-87: λ ≈ 698 nm = 6.98e-7 m
        let wl = FrequencyStandard::StrontiumOptical.transition_wavelength_m();
        assert!((wl - 6.98e-7).abs() < 1e-8);
    }

    #[test]
    fn ytterbium_wavelength() {
        // Yb-171: λ ≈ 578 nm = 5.78e-7 m
        let wl = FrequencyStandard::YtterbiumOptical.transition_wavelength_m();
        assert!((wl - 5.78e-7).abs() < 1e-8);
    }

    #[test]
    fn hydrogen_maser_wavelength() {
        // H maser: λ = c / 1.42e9 ≈ 21 cm
        let wl = FrequencyStandard::HydrogenMaser.transition_wavelength_m();
        assert!((wl - 0.211).abs() < 0.005);
    }

    #[test]
    fn quality_factors_reasonable() {
        assert!(FrequencyStandard::Cesium133.quality_factor() > 1e9);
        assert!(FrequencyStandard::StrontiumOptical.quality_factor() > 1e16);
        assert!(
            FrequencyStandard::HydrogenMaser.quality_factor()
                < FrequencyStandard::Cesium133.quality_factor()
        );
    }

    #[test]
    fn fractional_stability_ordering() {
        // Optical clocks should be more stable than microwave
        assert!(
            FrequencyStandard::StrontiumOptical.fractional_stability()
                < FrequencyStandard::Cesium133.fractional_stability()
        );
        assert!(
            FrequencyStandard::YtterbiumOptical.fractional_stability()
                < FrequencyStandard::Rubidium87.fractional_stability()
        );
    }

    #[test]
    fn all_standards_count() {
        assert_eq!(all_frequency_standards().len(), 5);
    }

    #[test]
    fn all_standards_unique_frequencies() {
        let stds = all_frequency_standards();
        for i in 0..stds.len() {
            for j in (i + 1)..stds.len() {
                assert!(
                    (stds[i].transition_frequency_hz() - stds[j].transition_frequency_hz()).abs()
                        > 1.0
                );
            }
        }
    }

    #[test]
    fn serde_roundtrip_frequency_standard() {
        for &std in all_frequency_standards() {
            let json = serde_json::to_string(&std).unwrap();
            let back: FrequencyStandard = serde_json::from_str(&json).unwrap();
            assert_eq!(std, back);
        }
    }

    // ── TimeScale ───────────────────────────────────────────────────────

    #[test]
    fn serde_roundtrip_time_scale() {
        let scales = [
            TimeScale::TAI,
            TimeScale::UTC,
            TimeScale::TT,
            TimeScale::GPS,
            TimeScale::TCB,
            TimeScale::TCG,
        ];
        for ts in &scales {
            let json = serde_json::to_string(ts).unwrap();
            let back: TimeScale = serde_json::from_str(&json).unwrap();
            assert_eq!(*ts, back);
        }
    }

    // ── Time scale conversions ──────────────────────────────────────────

    #[test]
    fn tai_tt_roundtrip() {
        let tai = 1_000_000.0;
        let tt = tai_to_tt(tai);
        let back = tt_to_tai(tt);
        assert!((back - tai).abs() < 1e-12);
    }

    #[test]
    fn tai_tt_offset() {
        assert!((tai_to_tt(0.0) - 32.184).abs() < f64::EPSILON);
    }

    #[test]
    fn tai_gps_roundtrip() {
        let tai = 500_000.0;
        let gps = tai_to_gps(tai);
        let back = gps_to_tai(gps);
        assert!((back - tai).abs() < 1e-12);
    }

    #[test]
    fn tai_gps_offset() {
        assert!((tai_to_gps(19.0)).abs() < f64::EPSILON);
    }

    #[test]
    fn tt_to_gps_via_tai() {
        // TT = TAI + 32.184, GPS = TAI - 19
        // So TT - GPS = 51.184
        let tai = 1_000_000.0;
        let tt = tai_to_tt(tai);
        let gps = tai_to_gps(tai);
        assert!((tt - gps - 51.184).abs() < 1e-10);
    }

    // ── Leap seconds ────────────────────────────────────────────────────

    #[test]
    fn leap_seconds_before_1972() {
        assert_eq!(leap_seconds_at(1971, 12), 0);
        assert_eq!(leap_seconds_at(1960, 6), 0);
    }

    #[test]
    fn leap_seconds_1972_jan() {
        assert_eq!(leap_seconds_at(1972, 1), 10);
    }

    #[test]
    fn leap_seconds_1972_jul() {
        assert_eq!(leap_seconds_at(1972, 7), 11);
    }

    #[test]
    fn leap_seconds_1972_between() {
        // Between Jan and Jul 1972, still 10
        assert_eq!(leap_seconds_at(1972, 6), 10);
    }

    #[test]
    fn leap_seconds_2017() {
        assert_eq!(leap_seconds_at(2017, 1), 37);
    }

    #[test]
    fn leap_seconds_2024() {
        // No leap seconds after 2017-01-01; ΔAT remains 37
        assert_eq!(leap_seconds_at(2024, 6), 37);
    }

    #[test]
    fn leap_seconds_gps_epoch() {
        // GPS epoch: 1980-01-06, ΔAT = 19
        assert_eq!(leap_seconds_at(1980, 1), 19);
    }

    #[test]
    fn tai_to_utc_offset_2024() {
        assert!((tai_to_utc_offset(2024, 1) - 37.0).abs() < f64::EPSILON);
    }

    #[test]
    fn utc_to_tai_offset_2024() {
        assert!((utc_to_tai_offset(2024, 1) - (-37.0)).abs() < f64::EPSILON);
    }

    #[test]
    fn tai_utc_offsets_inverse() {
        let year = 2020;
        let month = 6;
        assert!(
            (tai_to_utc_offset(year, month) + utc_to_tai_offset(year, month)).abs() < f64::EPSILON
        );
    }

    // ── AtomicInstant ───────────────────────────────────────────────────

    #[test]
    fn atomic_instant_basic() {
        let t = AtomicInstant::new(100, 500_000_000);
        assert!((t.tai_seconds() - 100.5).abs() < 1e-9);
    }

    #[test]
    fn atomic_instant_nanos_overflow() {
        // 1_500_000_000 ns = 1 s + 500_000_000 ns
        let t = AtomicInstant::new(10, 1_500_000_000);
        assert_eq!(t.seconds, 11);
        assert_eq!(t.nanos, 500_000_000);
    }

    #[test]
    fn atomic_instant_from_tai_seconds() {
        let t = AtomicInstant::from_tai_seconds(42.75);
        assert_eq!(t.seconds, 42);
        assert_eq!(t.nanos, 750_000_000);
    }

    #[test]
    fn atomic_instant_tai_seconds_roundtrip() {
        let s = 123_456_789.123_456_79;
        let t = AtomicInstant::from_tai_seconds(s);
        let back = t.tai_seconds();
        // Precision limited by f64 → i64 + u32 split
        assert!((back - s).abs() < 1e-6);
    }

    #[test]
    fn atomic_instant_to_tt() {
        let t = AtomicInstant::new(1000, 0);
        assert!((t.to_tt_seconds() - 1032.184).abs() < 1e-10);
    }

    #[test]
    fn atomic_instant_to_gps() {
        let t = AtomicInstant::new(1000, 0);
        assert!((t.to_gps_seconds() - 981.0).abs() < 1e-10);
    }

    #[test]
    fn atomic_instant_duration_since() {
        let t1 = AtomicInstant::new(100, 250_000_000);
        let t2 = AtomicInstant::new(200, 750_000_000);
        let dt = t2.duration_since(&t1);
        assert!((dt - 100.5).abs() < 1e-9);
    }

    #[test]
    fn atomic_instant_duration_since_negative() {
        let t1 = AtomicInstant::new(200, 0);
        let t2 = AtomicInstant::new(100, 0);
        let dt = t2.duration_since(&t1);
        assert!((dt - (-100.0)).abs() < 1e-9);
    }

    #[test]
    fn atomic_instant_add_seconds() {
        let t = AtomicInstant::new(100, 0);
        let t2 = t.add_seconds(1.5);
        assert!((t2.tai_seconds() - 101.5).abs() < 1e-9);
    }

    #[test]
    fn atomic_instant_add_negative() {
        let t = AtomicInstant::new(100, 0);
        let t2 = t.add_seconds(-10.0);
        assert!((t2.tai_seconds() - 90.0).abs() < 1e-9);
    }

    #[test]
    fn atomic_instant_ordering() {
        let t1 = AtomicInstant::new(100, 0);
        let t2 = AtomicInstant::new(100, 1);
        let t3 = AtomicInstant::new(101, 0);
        assert!(t1 < t2);
        assert!(t2 < t3);
    }

    #[test]
    fn serde_roundtrip_atomic_instant() {
        let t = AtomicInstant::new(1_700_000_000, 123_456_789);
        let json = serde_json::to_string(&t).unwrap();
        let back: AtomicInstant = serde_json::from_str(&json).unwrap();
        assert_eq!(t, back);
    }

    // ── Relativistic clock corrections ──────────────────────────────────

    #[test]
    fn gravitational_redshift_sign() {
        // Higher receiver: negative shift (blueshift)
        let shift = gravitational_redshift(100.0);
        assert!(shift < 0.0);
    }

    #[test]
    fn gravitational_redshift_zero() {
        assert!((gravitational_redshift(0.0)).abs() < f64::EPSILON);
    }

    #[test]
    fn gravitational_redshift_magnitude() {
        // 1 meter height: Δf/f ≈ -1.09e-16
        let shift = gravitational_redshift(1.0);
        assert!((shift - (-1.09e-16)).abs() < 1e-18);
    }

    #[test]
    fn gps_clock_correction() {
        // GPS: R ≈ 26,560 km, v ≈ 3874 m/s
        let correction = schwarzschild_clock_correction_us_per_day(26_560_000.0, 3874.0);
        // Expected: ~+38.6 μs/day
        assert!(
            (correction - 38.6).abs() < 1.0,
            "GPS correction = {correction} μs/day, expected ~38.6"
        );
    }

    #[test]
    fn gps_gravitational_part() {
        // Gravitational: GM/c² × (1/R_earth - 1/R_orbit) × 86400 × 1e6
        let c2 = C * C;
        let grav = GM_EARTH / c2 * (1.0 / 6_371_000.0 - 1.0 / 26_560_000.0) * 86400.0 * 1e6;
        assert!(
            (grav - 45.7).abs() < 0.5,
            "gravitational part = {grav} μs/day, expected ~45.7"
        );
    }

    #[test]
    fn gps_velocity_part() {
        // Velocity only: -v²/(2c²) × 86400 × 1e6
        let c2 = C * C;
        let vel = 3874.0 * 3874.0 / (2.0 * c2) * 86400.0 * 1e6;
        assert!(
            (vel - 7.2).abs() < 0.1,
            "velocity part = {vel} μs/day, expected ~7.2"
        );
    }

    #[test]
    fn schwarzschild_surface_zero_velocity() {
        // At Earth surface with zero velocity → no difference from ground clock
        let earth_radius = 6.371e6;
        let correction = schwarzschild_clock_correction_us_per_day(earth_radius, 0.0);
        assert!(correction.abs() < 0.01);
    }

    #[test]
    fn schwarzschild_leo() {
        // Low Earth orbit: R ≈ 6771 km (400 km altitude), v ≈ 7660 m/s
        // Velocity effect dominates → net correction should be negative
        let correction = schwarzschild_clock_correction_us_per_day(6_771_000.0, 7660.0);
        assert!(
            correction < 0.0,
            "LEO correction = {correction} μs/day, expected negative"
        );
    }

    #[test]
    fn second_order_doppler_zero() {
        assert!((second_order_doppler_shift(0.0)).abs() < f64::EPSILON);
    }

    #[test]
    fn second_order_doppler_sign() {
        // Moving clock runs slow → negative shift
        assert!(second_order_doppler_shift(0.1) < 0.0);
    }

    #[test]
    fn second_order_doppler_magnitude() {
        // β = 0.01 → Δf/f = -5e-5
        let shift = second_order_doppler_shift(0.01);
        assert!((shift - (-5e-5)).abs() < 1e-10);
    }

    #[test]
    fn sagnac_equator() {
        // At equator (latitude = 0), Sagnac correction = 0
        let dt = sagnac_correction_ns(0.0, 1e6);
        assert!(dt.abs() < 1e-10);
    }

    #[test]
    fn sagnac_pole() {
        // At north pole (latitude = π/2), maximum effect
        let dt = sagnac_correction_ns(core::f64::consts::FRAC_PI_2, 1e6);
        assert!(dt > 0.0);
        // 4 × 7.29e-5 × 1e6 × 1 / c² × 1e9
        let expected = 4.0 * EARTH_ROTATION_RAD_S * 1e6 / (C * C) * 1e9;
        assert!((dt - expected).abs() < 1e-6);
    }

    #[test]
    fn sagnac_sign_hemisphere() {
        // Northern hemisphere: positive
        let dt_north = sagnac_correction_ns(0.5, 1e6);
        assert!(dt_north > 0.0);
        // Southern hemisphere: negative
        let dt_south = sagnac_correction_ns(-0.5, 1e6);
        assert!(dt_south < 0.0);
    }

    // ── Constants ───────────────────────────────────────────────────────

    #[test]
    fn lg_rate_value() {
        assert!((LG_RATE - 6.969_290_134e-10).abs() < 1e-19);
    }

    #[test]
    fn lc_rate_value() {
        assert!((LC_RATE - 1.480_826_867_41e-8).abs() < 1e-18);
    }

    #[test]
    fn tai_tt_offset_value() {
        assert!((TAI_TT_OFFSET_S - 32.184).abs() < f64::EPSILON);
    }

    #[test]
    fn tai_gps_offset_value() {
        assert!((TAI_GPS_OFFSET_S - 19.0).abs() < f64::EPSILON);
    }

    // ── Leap second table completeness ──────────────────────────────────

    #[test]
    fn leap_second_table_length() {
        assert_eq!(LEAP_SECONDS.len(), 28);
    }

    #[test]
    fn leap_seconds_monotonic() {
        for i in 1..LEAP_SECONDS.len() {
            let (y_prev, m_prev, d_prev) = LEAP_SECONDS[i - 1];
            let (y_curr, m_curr, d_curr) = LEAP_SECONDS[i];
            assert!(
                (y_curr, m_curr) > (y_prev, m_prev),
                "Table not sorted at index {i}"
            );
            assert!(
                d_curr > d_prev,
                "ΔAT not monotonically increasing at index {i}: {d_prev} -> {d_curr}"
            );
        }
    }

    #[test]
    fn leap_seconds_first_and_last() {
        assert_eq!(LEAP_SECONDS[0], (1972, 1, 10));
        assert_eq!(LEAP_SECONDS[LEAP_SECONDS.len() - 1], (2017, 1, 37));
    }

    #[test]
    fn leap_seconds_all_valid_months() {
        for &(_, m, _) in LEAP_SECONDS {
            assert!(m == 1 || m == 7, "Leap second in unexpected month: {m}");
        }
    }

    // ── Cross-validation ────────────────────────────────────────────────

    #[test]
    fn cesium_defines_si_second() {
        // The SI second is exactly 9_192_631_770 periods of Cs-133
        let f = FrequencyStandard::Cesium133.transition_frequency_hz();
        assert!((f - 9_192_631_770.0).abs() < f64::EPSILON);
    }

    #[test]
    fn gps_epoch_leap_seconds_consistent() {
        // GPS epoch (1980-01-06): TAI-UTC = 19 s = TAI_GPS_OFFSET_S
        let delta_at = leap_seconds_at(1980, 1);
        assert!((delta_at as f64 - TAI_GPS_OFFSET_S).abs() < f64::EPSILON);
    }

    #[test]
    fn atomic_instant_tt_gps_consistent() {
        let t = AtomicInstant::new(1_000_000, 0);
        let tt = t.to_tt_seconds();
        let gps = t.to_gps_seconds();
        // TT - GPS = 32.184 + 19 = 51.184
        assert!((tt - gps - 51.184).abs() < 1e-10);
    }
}
