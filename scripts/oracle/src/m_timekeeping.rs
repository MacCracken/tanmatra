//! `timekeeping.tsv` — the port plan M5 grid for `tanmatra::timekeeping`.
//!
//! - **Leap lookups** for years 1955–2035 and months 0–13.
//! - **TAI − UTC by MJD:** a 3.7-day sweep from 1960 to 2026, every integer MJD where the value
//!   changes and 1e-6 before it, and NaN / −∞. (+∞ and huge MJDs overflow the oracle's integer
//!   calendar arithmetic, which panics under the oracle's overflow checks.)
//! - **AtomicInstant** at epoch 0, ±2.1e9 s and the ±9.2e18 s and i64 MIN / MAX edges, plus
//!   **exact .5 ns rounding ties** in `from_tai_seconds` and `add_seconds`. The ties are found by
//!   search, so that `(s − floor(s))·1e9` or `(dt − trunc(dt))·1e9` is exactly n + ½.
//! - **J2000 and T0** in the TCG / TCB conversions; the frequency standards; the relativistic
//!   corrections.
//!
//! An instant is the two tokens `i:seconds i:nanos`. Its private fields are read back exactly
//! through the public `nanoseconds_since(epoch)`. `nanoseconds_since` itself (an i128) is
//! `i:hi i:lo`: the high 64 bits signed, then the low 64 bits as an i64 bit pattern.
//! Classes: integer and table results are A, the arithmetic B, and `sagnac_correction_ns`
//! (sin) C.

use crate::tsv::{Class, Module, Tok, e, f, i};
use tanmatra::timekeeping::{self as tk, AtomicInstant, FrequencyStandard, TimeScale};

pub const STANDARDS: [FrequencyStandard; 5] = [
    FrequencyStandard::Cesium133,
    FrequencyStandard::Rubidium87,
    FrequencyStandard::HydrogenMaser,
    FrequencyStandard::StrontiumOptical,
    FrequencyStandard::YtterbiumOptical,
];

pub const SCALES: [TimeScale; 6] = [
    TimeScale::TAI,
    TimeScale::UTC,
    TimeScale::TT,
    TimeScale::GPS,
    TimeScale::TCB,
    TimeScale::TCG,
];

const NS: i128 = 1_000_000_000;

/// (seconds, nanos) of an instant, through the public API.
fn parts(t: &AtomicInstant) -> (i64, i64) {
    let v = t.nanoseconds_since(&AtomicInstant::new(0, 0));
    let s = i64::try_from(v.div_euclid(NS)).unwrap_or(0);
    let n = i64::try_from(v.rem_euclid(NS)).unwrap_or(0);
    (s, n)
}

fn inst(t: &AtomicInstant) -> [Tok; 2] {
    let (s, n) = parts(t);
    [i(s), i(n)]
}

/// The i128 as (signed high 64 bits, low 64 bits as an i64 bit pattern).
#[allow(clippy::cast_possible_truncation)]
fn i128_toks(v: i128) -> [Tok; 2] {
    [i((v >> 64) as i64), i(v as i64)]
}

/// Fractions whose ×1e9 is exactly n + ½, offset by each whole part (ties for round()).
fn tie_seconds() -> Vec<f64> {
    let mut v = Vec::new();
    for whole in [0.0, 3.0, 1000.0, 1e6, -1.0, -42.0] {
        let mut found = 0;
        let mut n = 0u32;
        while found < 6 && n < 2_000_000 {
            let s: f64 = whole + (f64::from(n) + 0.5) / 1e9;
            let frac = (s - libm::floor(s)) * 1e9;
            if frac - libm::floor(frac) == 0.5 {
                v.push(s);
                found += 1;
            }
            n += 7919;
        }
    }
    v
}

/// As `tie_seconds`, for add_seconds' truncation: (dt − trunc(dt))·1e9 = ±(n + ½).
fn tie_offsets() -> Vec<f64> {
    let mut v = Vec::new();
    for whole in [0.0, 2.0, -5.0, 1e5] {
        for sign in [1.0, -1.0] {
            let mut found = 0;
            let mut n = 0u32;
            while found < 4 && n < 2_000_000 {
                let dt: f64 = sign * (whole + (f64::from(n) + 0.5) / 1e9);
                let frac = libm::fabs((dt - libm::trunc(dt)) * 1e9);
                if frac - libm::floor(frac) == 0.5 {
                    v.push(dt);
                    found += 1;
                }
                n += 7919;
            }
        }
    }
    v
}

pub fn emit() -> Module {
    let mut m = Module::new("timekeeping");

    // ── constants ──
    for (name, v) in [
        ("TAI_TT_OFFSET_S", tk::TAI_TT_OFFSET_S),
        ("TAI_GPS_OFFSET_S", tk::TAI_GPS_OFFSET_S),
        ("LG_RATE", tk::LG_RATE),
        ("LC_RATE", tk::LC_RATE),
        ("LB_RATE", tk::LB_RATE),
        ("TDB0_S", tk::TDB0_S),
        ("T0_JD", tk::T0_JD),
    ] {
        m.row("const", Class::A, &[crate::tsv::s(name)], &[f(v)]);
    }
    let (vy, vm, vd) = tk::LEAP_SECOND_TABLE_VALID_UNTIL;
    m.row(
        "const",
        Class::A,
        &[crate::tsv::s("LEAP_SECOND_TABLE_VALID_UNTIL")],
        &[i(vy), i(vm), i(vd)],
    );

    // ── frequency standards ──
    let mut all = vec![Tok::V(tk::all_frequency_standards().len())];
    all.extend(
        tk::all_frequency_standards()
            .iter()
            .map(|&s| e("FrequencyStandard", s)),
    );
    m.row("all_frequency_standards", Class::A, &[], &all);
    for s in STANDARDS {
        let v = [e("FrequencyStandard", s)];
        m.row(
            "FrequencyStandard::transition_frequency_hz",
            Class::A,
            &v,
            &[f(s.transition_frequency_hz())],
        );
        m.row(
            "FrequencyStandard::transition_wavelength_m",
            Class::B,
            &v,
            &[f(s.transition_wavelength_m())],
        );
        m.row(
            "FrequencyStandard::fractional_stability",
            Class::A,
            &v,
            &[f(s.fractional_stability())],
        );
        for lw in [1.0, 0.5, 1e-3, 7.0, 0.0, -1.0, f64::INFINITY, f64::NAN] {
            let q = s.quality_factor_for_linewidth(lw);
            m.row(
                "FrequencyStandard::quality_factor_for_linewidth",
                Class::B,
                &[e("FrequencyStandard", s), f(lw)],
                &[f(q)],
            );
        }
    }

    // ── time scales ──
    let jds = [
        0.0,
        2_400_000.5,
        tk::T0_JD,
        2_443_144.5,
        2_451_545.0,
        2_460_000.5,
        2_488_069.5,
        2_816_787.5,
        -1e6,
        1e7,
    ];
    for &jd in &jds {
        m.row(
            "tcg_minus_tt_seconds",
            Class::B,
            &[f(jd)],
            &[f(tk::tcg_minus_tt_seconds(jd))],
        );
        m.row(
            "tt_to_tcg_jd",
            Class::B,
            &[f(jd)],
            &[f(tk::tt_to_tcg_jd(jd))],
        );
        m.row(
            "tcg_to_tt_jd",
            Class::B,
            &[f(jd)],
            &[f(tk::tcg_to_tt_jd(jd))],
        );
        m.row(
            "tcb_to_tdb_jd",
            Class::B,
            &[f(jd)],
            &[f(tk::tcb_to_tdb_jd(jd))],
        );
        m.row(
            "tdb_to_tcb_jd",
            Class::B,
            &[f(jd)],
            &[f(tk::tdb_to_tcb_jd(jd))],
        );
        m.row(
            "tcb_minus_tcg_secular_seconds",
            Class::B,
            &[f(jd)],
            &[f(tk::tcb_minus_tcg_secular_seconds(jd))],
        );
    }
    for t in [
        0.0,
        1.0,
        -1.0,
        0.5,
        2.1e9,
        1e18,
        -1e18,
        -0.0,
        f64::INFINITY,
        f64::NAN,
    ] {
        m.row("tai_to_tt", Class::B, &[f(t)], &[f(tk::tai_to_tt(t))]);
        m.row("tt_to_tai", Class::B, &[f(t)], &[f(tk::tt_to_tai(t))]);
        m.row("tai_to_gps", Class::B, &[f(t)], &[f(tk::tai_to_gps(t))]);
        m.row("gps_to_tai", Class::B, &[f(t)], &[f(tk::gps_to_tai(t))]);
    }

    // ── leap seconds ──
    for year in 1955..=2035 {
        for month in 0..=13u8 {
            let ins = [i(year), i(month)];
            m.row(
                "leap_seconds_at",
                Class::A,
                &ins,
                &[i(tk::leap_seconds_at(year, month))],
            );
            m.row(
                "tai_to_utc_offset",
                Class::A,
                &ins,
                &[f(tk::tai_to_utc_offset(year, month))],
            );
            m.row(
                "utc_to_tai_offset",
                Class::A,
                &ins,
                &[f(tk::utc_to_tai_offset(year, month))],
            );
        }
    }
    let mut mjds: Vec<f64> = (0..=6600).map(|k| 36_934.0 + f64::from(k) * 3.7).collect();
    let mut prev = tk::tai_minus_utc_seconds_mjd(36_000.0);
    for d in 36_001..=61_500 {
        let x = f64::from(d);
        let v = tk::tai_minus_utc_seconds_mjd(x);
        if v.to_bits() != prev.to_bits() && (v - prev).abs() > 1e-3 {
            mjds.extend([x - 1e-6, x]);
        }
        prev = v;
    }
    mjds.extend([f64::NAN, f64::NEG_INFINITY, 0.0, -1.0]);
    for &mjd in &mjds {
        m.row(
            "tai_minus_utc_seconds_mjd",
            Class::B,
            &[f(mjd)],
            &[f(tk::tai_minus_utc_seconds_mjd(mjd))],
        );
    }

    // ── AtomicInstant ──
    let secs = [
        0i64,
        1,
        -1,
        2_100_000_000,
        -2_100_000_000,
        9_200_000_000_000_000_000,
        -9_200_000_000_000_000_000,
        i64::MAX,
        i64::MIN,
    ];
    for &s in &secs {
        for n in [0u32, 1, 999_999_999, 1_000_000_000, 1_500_000_000, u32::MAX] {
            m.row(
                "AtomicInstant::new",
                Class::A,
                &[i(s), i(n)],
                &inst(&AtomicInstant::new(s, n)),
            );
        }
    }
    let mut ss = tie_seconds();
    ss.extend([
        0.0,
        0.25,
        1.5,
        -1.5,
        -0.0,
        5e-324,
        2.1e9 + 0.123_456_789,
        1e15 + 0.5,
        9.3e18,
        -9.3e18,
        1e300,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ]);
    for &s in &ss {
        m.row(
            "AtomicInstant::from_tai_seconds",
            Class::A,
            &[f(s)],
            &inst(&AtomicInstant::from_tai_seconds(s)),
        );
    }
    let instants: Vec<AtomicInstant> = secs
        .iter()
        .flat_map(|&s| [0u32, 1, 500_000_000, 999_999_999].map(move |n| AtomicInstant::new(s, n)))
        .collect();
    for t in &instants {
        let ins = inst(t);
        m.row(
            "AtomicInstant::tai_seconds",
            Class::B,
            &ins,
            &[f(t.tai_seconds())],
        );
        m.row(
            "AtomicInstant::to_tt_seconds",
            Class::B,
            &ins,
            &[f(t.to_tt_seconds())],
        );
        m.row(
            "AtomicInstant::to_gps_seconds",
            Class::B,
            &ins,
            &[f(t.to_gps_seconds())],
        );
        for dt in [
            0i64,
            1,
            -1,
            999_999_999,
            1_000_000_000,
            -1_000_000_001,
            123_456_789_012,
            i64::MAX,
            i64::MIN,
        ] {
            let mut a = ins.to_vec();
            a.push(i(dt));
            m.row(
                "AtomicInstant::add_nanoseconds",
                Class::A,
                &a,
                &inst(&t.add_nanoseconds(dt)),
            );
        }
        let mut dts = tie_offsets();
        dts.extend([
            0.0,
            1.5,
            -2.25,
            1e-9,
            4e-10,
            0.3,
            -0.3,
            1e18,
            -1e18,
            9.5e18,
            -9.5e18,
            1e300,
            f64::INFINITY,
            f64::NAN,
        ]);
        for dt in dts {
            let mut a = ins.to_vec();
            a.push(f(dt));
            m.row(
                "AtomicInstant::add_seconds",
                Class::A,
                &a,
                &inst(&t.add_seconds(dt)),
            );
        }
    }
    for p in &instants {
        for q in instants.iter().step_by(3) {
            let mut a = inst(p).to_vec();
            a.extend(inst(q));
            m.row(
                "AtomicInstant::nanoseconds_since",
                Class::A,
                &a,
                &i128_toks(p.nanoseconds_since(q)),
            );
            // duration_since subtracts the i64 seconds, which overflows (and panics under the
            // oracle's checks) across the extremes: only pairs whose difference fits.
            let (ps, _) = parts(p);
            let (qs, _) = parts(q);
            if ps.checked_sub(qs).is_some() {
                m.row(
                    "AtomicInstant::duration_since",
                    Class::B,
                    &a,
                    &[f(p.duration_since(q))],
                );
            }
        }
    }

    // ── relativistic corrections ──
    for h in [
        0.0,
        1.0,
        -1.0,
        100.0,
        8848.0,
        20_200_000.0,
        -0.0,
        1e300,
        f64::NAN,
    ] {
        m.row(
            "gravitational_redshift",
            Class::B,
            &[f(h)],
            &[f(tk::gravitational_redshift(h))],
        );
    }
    let gm = tanmatra::constants::GM_EARTH;
    let gps_a = 26_561_750.0;
    let orbits = [
        (gps_a, libm::sqrt(gm / gps_a)),
        (6_771_000.0, 7_660.0),
        (42_164_000.0, 3_075.0),
        (6_378_137.0, 0.0),
        (1e9, 1.0),
        (0.0, 0.0),
    ];
    for (r, v) in orbits {
        m.row(
            "schwarzschild_clock_correction_us_per_day",
            Class::B,
            &[f(r), f(v)],
            &[f(tk::schwarzschild_clock_correction_us_per_day(r, v))],
        );
    }
    for beta in [
        0.0,
        1e-9,
        1e-4,
        0.01,
        0.1,
        0.5,
        0.9,
        0.999_999,
        1.0,
        -1.0,
        1.5,
        -0.3,
        f64::NAN,
    ] {
        m.row(
            "second_order_doppler_shift",
            Class::B,
            &[f(beta)],
            &[f(tk::second_order_doppler_shift(beta))],
        );
        m.row(
            "time_dilation_shift_exact",
            Class::B,
            &[f(beta)],
            &[f(tk::time_dilation_shift_exact(beta))],
        );
    }
    let pi = core::f64::consts::PI;
    for lat in [0.0, pi / 6.0, pi / 4.0, pi / 2.0, -pi / 3.0, 0.8, 1.0] {
        for area in [1.0, 100.0, 1e6] {
            m.row(
                "sagnac_correction_ns",
                Class::C,
                &[f(lat), f(area)],
                &[f(tk::sagnac_correction_ns(lat, area))],
            );
        }
    }
    let re = 6_378_137.0;
    let pts = [
        (re, 0.0, 0.0, re),
        (0.0, re, re, 0.0),
        (re, 0.0, -re, 0.0),
        (gps_a, 0.0, 0.0, gps_a),
        (1.0, 2.0, 3.0, 4.0),
        (0.0, 0.0, 0.0, 0.0),
    ];
    for (x1, y1, x2, y2) in pts {
        m.row(
            "sagnac_one_way_ns",
            Class::B,
            &[f(x1), f(y1), f(x2), f(y2)],
            &[f(tk::sagnac_one_way_ns(x1, y1, x2, y2))],
        );
    }
    m
}
