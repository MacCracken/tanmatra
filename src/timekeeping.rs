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
//! - **CIPM 2025** (CCTF Recommendation 24-2): recommended frequencies of
//!   secondary representations of the SI second, in force since 2026-03-27
//! - **CODATA 2022**: Fundamental constants
//! - **IERS Bulletin C** / USNO `tai-utc.dat`: TAI − UTC, including 1961–1971
//! - **IERS Conventions (2010)**, Chapters 1 and 10: L_G, L_B, L_C, TDB0, T0

use crate::constants::{C, EARTH_ROTATION_RAD_S, GM_EARTH, STANDARD_GRAVITY};
use serde::{Deserialize, Serialize};

// ── Frequency Standards ─────────────────────────────────────────────────────

/// Atomic frequency standards used to define or approximate the SI second.
///
/// Each variant represents a real atomic transition used in precision
/// timekeeping. The cesium-133 hyperfine transition defines the SI second
/// (CGPM 1967, reaffirmed 2019). The Rb and optical standards are secondary
/// representations of the second (CIPM 2025 recommended values).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum FrequencyStandard {
    /// Cesium-133 hyperfine transition (SI second definition, CGPM 1967/2019).
    ///
    /// Transition: ⁶S₁/₂ (F=3 → F=4), frequency = 9 192 631 770 Hz (exact).
    Cesium133,
    /// Rubidium-87 hyperfine transition.
    ///
    /// Transition: ⁵S₁/₂ (F=1 → F=2), frequency 6 834 682 610.904 312 9 Hz.
    /// Reference: CIPM 2025 recommended value (u_r = 3.4e-16).
    Rubidium87,
    /// Hydrogen maser (1420 MHz, 21 cm line).
    ///
    /// Transition: 1S₁/₂ (F=0 → F=1), frequency 1 420 405 751.7667(9) Hz.
    /// Reference: Hellwig et al., IEEE Trans. Instrum. Meas. IM-19, 200 (1970).
    HydrogenMaser,
    /// Strontium-87 optical lattice clock (secondary representation of the second).
    ///
    /// Transition: ¹S₀ → ³P₀, frequency 429 228 004 229 872.992 Hz.
    /// Reference: CIPM 2025 recommended value (u_r = 1.7e-16).
    StrontiumOptical,
    /// Ytterbium-171 optical lattice clock (secondary representation of the second).
    ///
    /// Transition: ¹S₀ → ³P₀, frequency 518 295 836 590 863.632 Hz.
    /// Reference: CIPM 2025 recommended value (u_r = 1.7e-16).
    YtterbiumOptical,
}

impl FrequencyStandard {
    /// Transition frequency in Hz.
    ///
    /// For Cs-133, this is exact by definition (SI second). Rb, Sr and Yb are
    /// CIPM 2025 recommended values; H is Hellwig et al. (1970). The optical
    /// values are stored at f64 resolution (0.0625 Hz, 1.5e-16 relative),
    /// comparable to their recommended uncertainties.
    #[must_use]
    #[inline]
    #[allow(clippy::excessive_precision)] // published CIPM digits, rounded by f64
    pub fn transition_frequency_hz(self) -> f64 {
        match self {
            Self::Cesium133 => 9_192_631_770.0,
            Self::Rubidium87 => 6_834_682_610.904_312_9,
            Self::HydrogenMaser => 1_420_405_751.766_7,
            Self::StrontiumOptical => 429_228_004_229_872.992,
            Self::YtterbiumOptical => 518_295_836_590_863.632,
        }
    }

    /// Transition wavelength in meters (λ = c / f).
    #[must_use]
    #[inline]
    pub fn transition_wavelength_m(self) -> f64 {
        C / self.transition_frequency_hz()
    }

    /// Order-of-magnitude line quality factor for each clock technology.
    ///
    /// These are unsourced labels (10¹⁰ Cs and Rb, 10⁹ H maser, 10¹⁷ optical).
    /// The line Q of a real clock is set by its interrogation linewidth; compute
    /// it with [`FrequencyStandard::quality_factor_for_linewidth`].
    #[must_use]
    #[inline]
    #[deprecated(
        since = "1.3.0",
        note = "unsourced order-of-magnitude values; use quality_factor_for_linewidth"
    )]
    pub fn quality_factor(self) -> f64 {
        match self {
            Self::Cesium133 => 1e10,
            Self::Rubidium87 => 1e10,
            Self::HydrogenMaser => 1e9,
            Self::StrontiumOptical => 1e17,
            Self::YtterbiumOptical => 1e17,
        }
    }

    /// Line quality factor Q = ν / Δν for a given observed linewidth (FWHM, Hz).
    ///
    /// For Ramsey interrogation with free-evolution time T the fringe width is
    /// Δν ≈ 1/(2T). Returns 0.0 for a non-positive linewidth.
    #[must_use]
    #[inline]
    pub fn quality_factor_for_linewidth(self, linewidth_hz: f64) -> f64 {
        if linewidth_hz <= 0.0 {
            return 0.0;
        }
        self.transition_frequency_hz() / linewidth_hz
    }

    /// Representative fractional frequency instability (Allan deviation) at
    /// τ = 1 s for each clock technology.
    ///
    /// - Cs: commercial beam standard, Microchip 5071A standard tube ≤ 1.2e-11
    ///   (high-performance tube ≤ 5e-12).
    /// - Rb: commercial standard, SRS FS725 < 2e-11.
    /// - H maser: active maser, Microchip MHM-2020, 1.5e-13.
    /// - Sr lattice clock: 4.8e-17, two independent clocks (Oelker et al.,
    ///   Nat. Photon. 13, 714 (2019)).
    /// - Yb lattice clock: 1.4e-16 (Zhu et al., arXiv:2606.10514 (2026)).
    #[must_use]
    #[inline]
    pub fn fractional_stability(self) -> f64 {
        match self {
            Self::Cesium133 => 1.2e-11,
            Self::Rubidium87 => 2e-11,
            Self::HydrogenMaser => 1.5e-13,
            Self::StrontiumOptical => 4.8e-17,
            Self::YtterbiumOptical => 1.4e-16,
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

/// L_G (IAU 2000 Resolution B1.9, defining constant): dTT/dTCG = 1 − L_G.
pub const LG_RATE: f64 = 6.969_290_134e-10;

/// L_C (IERS Conventions 2010, Table 1.1): average rate of TCB relative to TCG,
/// ⟨dTCB/dTCG⟩ = 1 + L_C.
pub const LC_RATE: f64 = 1.480_826_867_41e-8;

/// L_B (IAU 2006 Resolution B3, defining constant): dTDB/dTCB = 1 − L_B.
pub const LB_RATE: f64 = 1.550_519_768e-8;

/// TDB0 (IAU 2006 Resolution B3, defining constant), in seconds.
pub const TDB0_S: f64 = -6.55e-5;

/// T0: 1977-01-01T00:00:00 TAI as a TT (= TCG = TCB) Julian date
/// (IERS Conventions 2010, eq. 10.1).
pub const T0_JD: f64 = 2_443_144.500_372_5;

/// TCG − TT in seconds for a TT Julian date (IERS Conventions 2010, eq. 10.1):
///
/// TCG − TT = (L_G / (1 − L_G)) × (JD_TT − T0) × 86400 s.
///
/// f64 Julian dates resolve ≈ 40 µs, far finer than this difference changes.
#[must_use]
#[inline]
pub fn tcg_minus_tt_seconds(jd_tt: f64) -> f64 {
    LG_RATE / (1.0 - LG_RATE) * (jd_tt - T0_JD) * 86_400.0
}

/// Converts a TT Julian date to a TCG Julian date (IERS Conventions 2010, eq. 10.1).
#[must_use]
#[inline]
pub fn tt_to_tcg_jd(jd_tt: f64) -> f64 {
    jd_tt + tcg_minus_tt_seconds(jd_tt) / 86_400.0
}

/// Converts a TCG Julian date to a TT Julian date (inverse of [`tt_to_tcg_jd`]):
/// JD_TT = JD_TCG − L_G (JD_TCG − T0).
#[must_use]
#[inline]
pub fn tcg_to_tt_jd(jd_tcg: f64) -> f64 {
    jd_tcg - LG_RATE * (jd_tcg - T0_JD)
}

/// Converts a TCB Julian date to a TDB Julian date (IAU 2006 Resolution B3,
/// IERS Conventions 2010, eq. 10.3):
///
/// TDB = TCB − L_B × (JD_TCB − T0) × 86400 s + TDB0.
#[must_use]
#[inline]
pub fn tcb_to_tdb_jd(jd_tcb: f64) -> f64 {
    jd_tcb - LB_RATE * (jd_tcb - T0_JD) + TDB0_S / 86_400.0
}

/// Converts a TDB Julian date to a TCB Julian date (exact inverse of
/// [`tcb_to_tdb_jd`]).
#[must_use]
#[inline]
pub fn tdb_to_tcb_jd(jd_tdb: f64) -> f64 {
    (jd_tdb - LB_RATE * T0_JD - TDB0_S / 86_400.0) / (1.0 - LB_RATE)
}

/// Secular part of TCB − TCG in seconds for a TT Julian date
/// (IERS Conventions 2010, eq. 10.5):
///
/// TCB − TCG ≈ L_C × (JD_TT − T0) × 86400 s / (1 − L_B).
///
/// Omitted: the periodic terms P(TT) − P(T0) (amplitude ≈ 1.6 ms, which need
/// a solar-system ephemeris) and the position-dependent term c⁻² v_e·(x − x_e).
#[must_use]
#[inline]
pub fn tcb_minus_tcg_secular_seconds(jd_tt: f64) -> f64 {
    LC_RATE * (jd_tt - T0_JD) * 86_400.0 / (1.0 - LB_RATE)
}

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
/// As of IERS Bulletin C 72 (2026-07-06), no further leap seconds have been
/// announced; TAI − UTC remains 37 s. The table is valid until
/// [`LEAP_SECOND_TABLE_VALID_UNTIL`].
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

/// Date (year, month, day) until which the leap-second table is known to be
/// complete: the expiry of the IANA `leap-seconds.list` published after
/// IERS Bulletin C 72 (2026-06-28 + 1 y → 2027-06-28).
pub const LEAP_SECOND_TABLE_VALID_UNTIL: (i32, u8, u8) = (2027, 6, 28);

/// TAI − UTC rubber-second segments 1961–1971 from the USNO `tai-utc.dat`
/// table: (MJD start, offset s, MJD reference, rate s/day), valid from the
/// start until the next segment. TAI − UTC = offset + (MJD − MJD_ref) × rate.
const TAI_UTC_PRE_1972: &[(f64, f64, f64, f64)] = &[
    (37_300.0, 1.422_818_0, 37_300.0, 0.001_296),
    (37_512.0, 1.372_818_0, 37_300.0, 0.001_296),
    (37_665.0, 1.845_858_0, 37_665.0, 0.001_123_2),
    (38_334.0, 1.945_858_0, 37_665.0, 0.001_123_2),
    (38_395.0, 3.240_130_0, 38_761.0, 0.001_296),
    (38_486.0, 3.340_130_0, 38_761.0, 0.001_296),
    (38_639.0, 3.440_130_0, 38_761.0, 0.001_296),
    (38_761.0, 3.540_130_0, 38_761.0, 0.001_296),
    (38_820.0, 3.640_130_0, 38_761.0, 0.001_296),
    (38_942.0, 3.740_130_0, 38_761.0, 0.001_296),
    (39_004.0, 3.840_130_0, 38_761.0, 0.001_296),
    (39_126.0, 4.313_170_0, 39_126.0, 0.002_592),
    (39_887.0, 4.213_170_0, 39_126.0, 0.002_592),
];

/// MJD of 1972-01-01, when integer leap seconds began (TAI − UTC = 10 s).
const MJD_1972: f64 = 41_317.0;

/// Returns TAI − UTC in seconds for a UTC modified Julian date, including the
/// fractional "rubber second" era 1961-01-01 to 1971-12-31 (USNO `tai-utc.dat`).
///
/// Returns 0.0 before 1961-01-01 (MJD 37300), when UTC was not defined
/// against TAI.
#[must_use]
pub fn tai_minus_utc_seconds_mjd(mjd_utc: f64) -> f64 {
    if mjd_utc < MJD_1972 {
        let mut value = 0.0;
        for &(start, offset, reference, rate) in TAI_UTC_PRE_1972 {
            if mjd_utc >= start {
                value = offset + (mjd_utc - reference) * rate;
            } else {
                break;
            }
        }
        return value;
    }
    // Integer era: convert MJD to a calendar (year, month).
    let (year, month) = mjd_to_year_month(mjd_utc);
    f64::from(leap_seconds_at(year, month))
}

/// Converts a modified Julian date to its Gregorian (year, month).
fn mjd_to_year_month(mjd: f64) -> (i32, u8) {
    // Fliegel & Van Flandern (1968) algorithm on the integer Julian day number.
    #[allow(clippy::cast_possible_truncation)]
    let jdn = libm::floor(mjd + 2_400_001.0) as i64;
    let l = jdn + 68_569;
    let n = 4 * l / 146_097;
    let l = l - (146_097 * n + 3) / 4;
    let i = 4000 * (l + 1) / 1_461_001;
    let l = l - 1461 * i / 4 + 31;
    let j = 80 * l / 2447;
    let l2 = j / 11;
    let month = j + 2 - 12 * l2;
    let year = 100 * (n - 49) + i + l2;
    #[allow(clippy::cast_sign_loss)]
    (year as i32, month as u8)
}

/// Returns TAI − UTC (ΔAT) in whole seconds for a given date in the
/// integer-leap-second era (from 1972-01-01).
///
/// Scans the leap second table to find the most recent entry at or before
/// the given (year, month). Returns 0 for dates before 1972-01-01, when
/// TAI − UTC was not an integer; use [`tai_minus_utc_seconds_mjd`] for
/// 1961–1971.
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
#[serde(from = "AtomicInstantRepr")]
pub struct AtomicInstant {
    /// Whole seconds since 1958-01-01T00:00:00 TAI.
    seconds: i64,
    /// Sub-second nanoseconds \[0, 999_999_999\].
    nanos: u32,
}

/// Serialized form of [`AtomicInstant`]; deserialization normalizes nanoseconds.
#[derive(Deserialize)]
struct AtomicInstantRepr {
    seconds: i64,
    nanos: u32,
}

impl From<AtomicInstantRepr> for AtomicInstant {
    fn from(raw: AtomicInstantRepr) -> Self {
        Self::new(raw.seconds, raw.nanos)
    }
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

    /// Creates an `AtomicInstant` from floating-point TAI seconds since epoch,
    /// rounded to the nearest nanosecond.
    ///
    /// An f64 near today's epoch (≈ 2e9 s) resolves only ≈ 0.5 µs, so this
    /// conversion cannot carry nanosecond information; use
    /// [`AtomicInstant::new`] or [`AtomicInstant::add_nanoseconds`] for that.
    #[must_use]
    pub fn from_tai_seconds(s: f64) -> Self {
        if !s.is_finite() {
            return Self::new(0, 0);
        }
        let whole = libm::floor(s);
        let frac_ns = libm::round((s - whole) * NANOS_PER_SEC);
        #[allow(clippy::cast_sign_loss)]
        let nanos = frac_ns.clamp(0.0, 1_000_000_000.0) as u32;
        Self::new(whole as i64, nanos)
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

    /// Returns the exact duration in nanoseconds between this instant and an
    /// earlier one (negative if `earlier` is later).
    #[must_use]
    pub fn nanoseconds_since(&self, earlier: &Self) -> i128 {
        (i128::from(self.seconds) - i128::from(earlier.seconds)) * 1_000_000_000
            + (i128::from(self.nanos) - i128::from(earlier.nanos))
    }

    /// Returns a new instant advanced by an exact number of nanoseconds
    /// (saturating at the i64 seconds range).
    #[must_use]
    pub fn add_nanoseconds(&self, dt_ns: i64) -> Self {
        let total = i128::from(self.nanos) + i128::from(dt_ns);
        let carry = total.div_euclid(1_000_000_000);
        let nanos = total.rem_euclid(1_000_000_000);
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let seconds = i128::from(self.seconds)
            .saturating_add(carry)
            .clamp(i128::from(i64::MIN), i128::from(i64::MAX)) as i64;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        Self {
            seconds,
            nanos: nanos as u32,
        }
    }

    /// Returns a new instant advanced by `dt` seconds.
    ///
    /// The offset is rounded to the nearest nanosecond and added in integer
    /// arithmetic, so the stored nanoseconds are never degraded by the size of
    /// the epoch. Non-finite `dt` returns `self` unchanged.
    #[must_use]
    pub fn add_seconds(&self, dt: f64) -> Self {
        if !dt.is_finite() {
            return *self;
        }
        let whole = libm::trunc(dt);
        let frac_ns = libm::round((dt - whole) * NANOS_PER_SEC);
        #[allow(clippy::cast_possible_truncation)]
        let whole_i = whole.clamp(-9.0e18, 9.0e18) as i64;
        #[allow(clippy::cast_possible_truncation)]
        let frac_i = frac_ns as i64;
        let shifted = Self {
            seconds: self.seconds.saturating_add(whole_i),
            nanos: self.nanos,
        };
        shifted.add_nanoseconds(frac_i)
    }
}

// ── Relativistic Clock Corrections ──────────────────────────────────────────

/// Gravitational redshift fractional frequency shift.
///
/// Δf/f = −g Δh / c², where g = 9.80665 m/s² (standard gravity).
///
/// A positive `delta_h_m` (receiver higher than emitter) gives a *negative*
/// fractional shift: light climbing out of the potential well is redshifted.
/// Equivalently, a clock raised by Δh runs *fast* by +gΔh/c² relative to the
/// lower clock. First order in Δh (uniform field g = 9.80665 m/s²).
///
/// # Arguments
///
/// * `delta_h_m` — Height difference in meters (positive = upward).
#[must_use]
#[inline]
pub fn gravitational_redshift(delta_h_m: f64) -> f64 {
    -STANDARD_GRAVITY * delta_h_m / (C * C)
}

/// Total relativistic clock-rate correction for a satellite in orbit relative
/// to a clock on the geoid (i.e. to TT).
///
/// A clock on the geoid runs at dτ/dTCG = 1 − L_G with L_G = W₀/c²
/// (W₀ = 62 636 856 m²/s², the geoid potential including Earth rotation). A
/// satellite clock runs at dτ/dTCG = 1 − (GM/r + v²/2)/c². The difference, in
/// µs/day, is
///
/// - Gravitational: (L_G − GM/(c² r)) × 86400 × 10⁶
/// - Velocity: −v²/(2c²) × 86400 × 10⁶
///
/// (IERS Conventions 2010, eq. 10.9; Ashby, Living Rev. Relativ. 6, 1 (2003)).
/// GM_earth = 3.986004418 × 10¹⁴ m³/s² (IERS 2010).
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
/// a = 26 561.75 km, v = √(GM/a) = 3873.8 m/s → gravitational +45.788 μs/day,
/// velocity −7.213 μs/day, net +38.575 μs/day (4.4647e-10).
#[must_use]
pub fn schwarzschild_clock_correction_us_per_day(
    orbital_radius_m: f64,
    orbital_velocity_m_s: f64,
) -> f64 {
    let c2 = C * C;
    let seconds_per_day: f64 = 86400.0;
    let us_per_s: f64 = 1e6;

    let grav = LG_RATE - GM_EARTH / (c2 * orbital_radius_m);
    let vel = orbital_velocity_m_s * orbital_velocity_m_s / (2.0 * c2);

    (grav - vel) * seconds_per_day * us_per_s
}

/// Second-order (transverse) Doppler shift, leading order in β.
///
/// Δf/f = −β²/2, where β = v/c.
///
/// This is the purely relativistic time dilation effect for a moving clock.
/// The exact value is √(1 − β²) − 1 ([`time_dilation_shift_exact`]); the
/// leading-order form is 2.5e-5 relatively off at β = 0.01 and 6.7% at β = 0.5.
///
/// # Arguments
///
/// * `beta` — v/c (dimensionless speed parameter).
#[must_use]
#[inline]
pub fn second_order_doppler_shift(beta: f64) -> f64 {
    -beta * beta / 2.0
}

/// Exact fractional frequency shift of a moving clock, √(1 − β²) − 1.
///
/// Computed as −β²/(1 + √(1 − β²)) to avoid cancellation at small β.
/// Returns −1.0 for |β| ≥ 1.
#[must_use]
#[inline]
pub fn time_dilation_shift_exact(beta: f64) -> f64 {
    let b2 = beta * beta;
    if b2.is_nan() || b2 >= 1.0 {
        return -1.0;
    }
    -b2 / (1.0 + libm::sqrt(1.0 - b2))
}

/// Sagnac time difference of a closed-loop interferometer (ring laser or
/// fibre gyroscope) lying horizontally on the rotating Earth, in nanoseconds.
///
/// Δt = 4 Ω A sin(φ) / c²
///
/// between the two counter-propagating beams, where A is the loop area and
/// Ω sin φ the component of Earth's rotation normal to it at latitude φ.
/// For one-way signal time transfer use [`sagnac_one_way_ns`].
///
/// # Arguments
///
/// * `latitude_rad` — Geodetic latitude in radians.
/// * `area_m2` — Area enclosed by the loop in m².
#[must_use]
pub fn sagnac_correction_ns(latitude_rad: f64, area_m2: f64) -> f64 {
    let c2 = C * C;
    let dt_s = 4.0 * EARTH_ROTATION_RAD_S * area_m2 * libm::sin(latitude_rad) / c2;
    dt_s * 1e9
}

/// Sagnac correction for a one-way signal between two points given in an
/// Earth-fixed frame, in nanoseconds.
///
/// Δt = 2 Ω A_z / c² = Ω (x₁ y₂ − y₁ x₂) / c²
///
/// where A_z is the area swept by the position vector projected on the
/// equatorial plane (Ashby, Living Rev. Relativ. 6, 1 (2003), eq. 32).
/// Positive when the signal travels eastward (with the rotation). A signal
/// once around the equator takes 207.4 ns longer eastward.
///
/// # Arguments
///
/// * `x1_m`, `y1_m` — equatorial components of the transmitter position (m).
/// * `x2_m`, `y2_m` — equatorial components of the receiver position (m).
#[must_use]
pub fn sagnac_one_way_ns(x1_m: f64, y1_m: f64, x2_m: f64, y2_m: f64) -> f64 {
    EARTH_ROTATION_RAD_S * (x1_m * y2_m - y1_m * x2_m) / (C * C) * 1e9
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
    #[allow(deprecated)]
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
        // Higher receiver: negative shift (redshift of light climbing out)
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
    fn schwarzschild_gps_iers_reference() {
        // IERS Conventions 2010 eq. 10.9 / Ashby 2003: a = 26 561.75 km,
        // net +38.575 µs/day (gravitational +45.788, velocity −7.213).
        let a = 26_561_750.0;
        let v = libm::sqrt(GM_EARTH / a);
        let net = schwarzschild_clock_correction_us_per_day(a, v);
        assert!((net - 38.575).abs() < 1e-3, "net={net}");
        let grav_only = schwarzschild_clock_correction_us_per_day(a, 0.0);
        assert!((grav_only - 45.788).abs() < 1e-3, "grav={grav_only}");
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
    fn sagnac_one_way_equatorial_loop() {
        // Eastward signal once around the equator: 2Ω(πR²)/c² = 207.39 ns,
        // summed from short one-way hops.
        let r = 6_378_137.0;
        let steps = 3600;
        let mut total = 0.0;
        for k in 0..steps {
            let a1 = 2.0 * core::f64::consts::PI * f64::from(k) / f64::from(steps);
            let a2 = 2.0 * core::f64::consts::PI * f64::from(k + 1) / f64::from(steps);
            total += sagnac_one_way_ns(
                r * libm::cos(a1),
                r * libm::sin(a1),
                r * libm::cos(a2),
                r * libm::sin(a2),
            );
        }
        let exact = 2.0 * EARTH_ROTATION_RAD_S * core::f64::consts::PI * r * r / (C * C) * 1e9;
        assert!((exact - 207.39).abs() < 0.01, "exact={exact}");
        assert!((total - exact).abs() / exact < 1e-5, "total={total}");
    }

    #[test]
    #[allow(clippy::excessive_precision)]
    fn cipm_2025_frequencies() {
        assert!(
            (FrequencyStandard::StrontiumOptical.transition_frequency_hz()
                - 429_228_004_229_872.992)
                .abs()
                < 0.07
        );
        assert!(
            (FrequencyStandard::YtterbiumOptical.transition_frequency_hz()
                - 518_295_836_590_863.632)
                .abs()
                < 0.07
        );
        assert!(
            (FrequencyStandard::Rubidium87.transition_frequency_hz() - 6_834_682_610.904_312_9)
                .abs()
                < 1e-6
        );
        let q = FrequencyStandard::Cesium133.quality_factor_for_linewidth(1.0);
        assert!((q - 9_192_631_770.0).abs() < 1e-3);
    }

    #[test]
    fn tcg_and_tdb_at_j2000() {
        let jd = 2_451_545.0;
        // TCG − TT at J2000 ≈ 0.505833 s.
        assert!((tcg_minus_tt_seconds(jd) - 0.505_833).abs() < 1e-6);
        assert!((tcg_to_tt_jd(tt_to_tcg_jd(jd)) - jd).abs() < 1e-9);
        // TDB − TCB at J2000 ≈ −11.2537 s.
        let tdb = tcb_to_tdb_jd(jd);
        assert!(
            ((tdb - jd) * 86_400.0 + 11.253_7).abs() < 1e-3,
            "{}",
            (tdb - jd) * 86_400.0
        );
        assert!((tdb_to_tcb_jd(tdb) - jd).abs() < 1e-9);
        assert!(tcg_minus_tt_seconds(T0_JD).abs() < 1e-12);
        assert!(tcb_minus_tcg_secular_seconds(T0_JD).abs() < 1e-12);
    }

    #[test]
    fn tai_minus_utc_rubber_second_era() {
        assert!((tai_minus_utc_seconds_mjd(37_300.0) - 1.422_818).abs() < 1e-9);
        // End of 1971: 4.2131700 + (41317 − 39126) × 0.002592 = 9.892242 s.
        assert!((tai_minus_utc_seconds_mjd(41_316.999) - 9.892_239).abs() < 1e-5);
        assert!((tai_minus_utc_seconds_mjd(41_317.0) - 10.0).abs() < 1e-12);
        assert!((tai_minus_utc_seconds_mjd(57_754.0) - 37.0).abs() < 1e-12); // 2017-01-01
        assert!((tai_minus_utc_seconds_mjd(57_753.0) - 36.0).abs() < 1e-12); // 2016-12-31
        assert!(tai_minus_utc_seconds_mjd(30_000.0).abs() < 1e-12);
    }

    #[test]
    fn atomic_instant_exact_at_modern_epoch() {
        let t0 = AtomicInstant::new(2_100_000_000, 123_456_789);
        assert_eq!(t0.add_seconds(0.0), t0);
        let t1 = t0.add_seconds(1e-9);
        assert_eq!(t1.nanoseconds_since(&t0), 1);
        let t2 = t0.add_seconds(-2.5);
        assert_eq!(t2.nanoseconds_since(&t0), -2_500_000_000);
        let t3 = t0.add_nanoseconds(900_000_000);
        assert_eq!((t3.seconds, t3.nanos), (2_100_000_001, 23_456_789));
        assert_eq!(AtomicInstant::from_tai_seconds(0.3).nanos, 300_000_000);
    }

    #[test]
    fn atomic_instant_serde_normalizes() {
        let t: AtomicInstant = serde_json::from_str(r#"{"seconds":5,"nanos":4000000000}"#).unwrap();
        assert_eq!((t.seconds, t.nanos), (9, 0));
        assert!(t > AtomicInstant::new(8, 0));
    }

    #[test]
    fn time_dilation_exact_vs_leading_order() {
        let exact = time_dilation_shift_exact(0.5);
        assert!((exact - (libm::sqrt(0.75) - 1.0)).abs() < 1e-15);
        assert!((second_order_doppler_shift(1e-4) - time_dilation_shift_exact(1e-4)).abs() < 1e-16);
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
