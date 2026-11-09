//! `src/timekeeping_data.cyr` — frequency standards, time-scale constants, the leap-second
//! table and the 1961–1971 TAI − UTC segments (port plan M5).
//!
//! - **The frequency standards and public constants:** through the oracle's API.
//! - **The leap-second table (private):** DERIVED from the public `leap_seconds_at`. Every
//!   change between consecutive months is one entry, and ΔAT only grows. It is proven by
//!   re-running `leap_seconds_at` from the derived table over years 1900–2200 and months 0–13.
//! - **The USNO segments and MJD_1972 (private):** extracted by `rust_src`, then proven by
//!   recomputing `tai_minus_utc_seconds_mjd` over MJD 36500–42500 in quarter days, the
//!   segment starts ±1e-6, and NaN / −∞.

use crate::cyr::{Elem, doc_comment, f64_lit, header, snake, table};
use crate::m_timekeeping::{SCALES, STANDARDS};
use crate::rust_src::{f64_const, f64_tuples};
use std::fmt::Write as _;
use tanmatra::timekeeping as tk;

/// (year, month, ΔAT after) from the public lookup.
fn derive_leap_table() -> Vec<(i32, u8, i32)> {
    let mut out = Vec::new();
    let mut prev = 0;
    for year in 1950..=2100 {
        for month in 1..=12u8 {
            let v = tk::leap_seconds_at(year, month);
            if v != prev {
                out.push((year, month, v));
                prev = v;
            }
        }
    }
    out
}

/// The oracle's `leap_seconds_at`, transcribed over a table.
fn leap_at(table: &[(i32, u8, i32)], year: i32, month: u8) -> i32 {
    let mut result = 0;
    for &(y, m, d) in table {
        if year > y || (year == y && month >= m) {
            result = d;
        } else {
            break;
        }
    }
    result
}

/// The oracle's `mjd_to_year_month` (Fliegel & Van Flandern), transcribed.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn mjd_ym(mjd: f64) -> (i32, u8) {
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
    (year as i32, month as u8)
}

/// The oracle's `tai_minus_utc_seconds_mjd`, transcribed (the proof of the segments).
fn tai_minus_utc(seg: &[Vec<f64>], mjd_1972: f64, table: &[(i32, u8, i32)], mjd: f64) -> f64 {
    if mjd < mjd_1972 {
        let mut value = 0.0;
        for s in seg {
            if mjd >= s[0] {
                value = s[1] + (mjd - s[2]) * s[3];
            } else {
                break;
            }
        }
        return value;
    }
    let (y, m) = mjd_ym(mjd);
    f64::from(leap_at(table, y, m))
}

fn member(e: &str, v: &str) -> String {
    format!("TM_{}_{}", snake(e), snake(v))
}

pub fn emit(src: &str, provenance: &str) -> Result<String, String> {
    // ── derive / extract, then prove ──
    let leaps = derive_leap_table();
    for year in 1900..=2200 {
        for month in 0..=13u8 {
            if leap_at(&leaps, year, month) != tk::leap_seconds_at(year, month) {
                return Err(format!("leap-table proof failed at ({year}, {month})"));
            }
        }
    }
    let seg = f64_tuples(src, "TAI_UTC_PRE_1972", 4)?;
    let mjd_1972 = f64_const(src, "MJD_1972")?;
    let mut probes: Vec<f64> = (0..=24_000)
        .map(|k| 36_500.0 + f64::from(k) * 0.25)
        .collect();
    for s in &seg {
        probes.extend([s[0] - 1e-6, s[0], s[0] + 1e-6]);
    }
    probes.extend([mjd_1972 - 1e-6, mjd_1972, f64::NAN, f64::NEG_INFINITY]);
    for &mjd in &probes {
        let got = tai_minus_utc(&seg, mjd_1972, &leaps, mjd).to_bits();
        let want = tk::tai_minus_utc_seconds_mjd(mjd).to_bits();
        if got != want {
            return Err(format!("USNO-segment proof failed at MJD {mjd}"));
        }
    }

    let mut out = header(
        "timekeeping data: frequency standards, time-scale constants, leap seconds, TAI − UTC 1961–1971",
        provenance,
        &[
            "The leap-second table is derived from the oracle's leap_seconds_at, and the 1961–1971 segments are\n\
           extracted from rust-old/src/timekeeping.rs. Both are proven against the oracle before this file is\n\
           written (scripts/oracle/src/t_timekeeping.rs). tests/golden/timekeeping.tsv pins every use.",
        ],
    );

    // ── enums ──
    for (rust, variants) in [
        (
            "FrequencyStandard",
            STANDARDS
                .iter()
                .map(|v| format!("{v:?}"))
                .collect::<Vec<_>>(),
        ),
        (
            "TimeScale",
            SCALES.iter().map(|v| format!("{v:?}")).collect::<Vec<_>>(),
        ),
    ] {
        out.push_str(&doc_comment(src, &format!("pub enum {rust} {{"))?);
        let _ = writeln!(out, "enum Tm{rust} {{");
        for (k, v) in variants.iter().enumerate() {
            let _ = writeln!(out, "    {} = {k};", member(rust, v));
        }
        out.push_str("}\n");
        let _ = writeln!(
            out,
            "const TM_{}_COUNT = {};\n",
            snake(rust),
            variants.len()
        );
    }

    // ── frequency standards ──
    let lab: Vec<String> = STANDARDS.iter().map(|v| format!("{v:?}")).collect();
    let f = STANDARDS
        .iter()
        .map(|s| s.transition_frequency_hz())
        .collect();
    out.push_str(&table(
        "TM_FREQUENCY_STANDARD_TRANSITION_FREQUENCY_HZ",
        "Transition frequency, Hz (Cs exact by definition; Rb, Sr, Yb CIPM 2025; H Hellwig et al. 1970).\n\
         Indexed by TmFrequencyStandard; mirrors `FrequencyStandard::transition_frequency_hz`.",
        &Elem::F64(f),
        &lab,
    )?);
    let st = STANDARDS.iter().map(|s| s.fractional_stability()).collect();
    out.push_str(&table(
        "TM_FREQUENCY_STANDARD_FRACTIONAL_STABILITY",
        "Fractional instability (Allan deviation) at τ = 1 s, sourced per technology.\n\
         Indexed by TmFrequencyStandard; mirrors `FrequencyStandard::fractional_stability`.",
        &Elem::F64(st),
        &lab,
    )?);

    // ── public constants ──
    for (name, v) in [
        ("TAI_TT_OFFSET_S", tk::TAI_TT_OFFSET_S),
        ("TAI_GPS_OFFSET_S", tk::TAI_GPS_OFFSET_S),
        ("LG_RATE", tk::LG_RATE),
        ("LC_RATE", tk::LC_RATE),
        ("LB_RATE", tk::LB_RATE),
        ("TDB0_S", tk::TDB0_S),
        ("T0_JD", tk::T0_JD),
    ] {
        out.push_str(&doc_comment(src, &format!("pub const {name}:"))?);
        let _ = writeln!(out, "const TM_{name} = {};\n", f64_lit(v)?);
    }
    let (vy, vm, vd) = tk::LEAP_SECOND_TABLE_VALID_UNTIL;
    out.push_str(&doc_comment(
        src,
        "pub const LEAP_SECOND_TABLE_VALID_UNTIL:",
    )?);
    let _ = writeln!(out, "const TM_LEAP_SECOND_TABLE_VALID_UNTIL_YEAR = {vy};");
    let _ = writeln!(out, "const TM_LEAP_SECOND_TABLE_VALID_UNTIL_MONTH = {vm};");
    let _ = writeln!(out, "const TM_LEAP_SECOND_TABLE_VALID_UNTIL_DAY = {vd};\n");

    // ── leap seconds ──
    let lab: Vec<String> = leaps
        .iter()
        .map(|(y, m, _)| format!("{y}-{m:02}"))
        .collect();
    out.push_str(&doc_comment(src, "const LEAP_SECONDS:")?);
    let _ = writeln!(out, "const _TM_LEAP_SECOND_COUNT = {};\n", leaps.len());
    let _ = writeln!(
        out,
        "# Column: year.\nvar _TM_LEAP_SECOND_YEAR: u16[{}] = {{",
        leaps.len()
    );
    for (y, _, _) in &leaps {
        let _ = writeln!(out, "    {y},");
    }
    out.push_str("};\n\n");
    let months = leaps.iter().map(|l| l.1).collect();
    out.push_str(&table(
        "_TM_LEAP_SECOND_MONTH",
        "Column: month (1 = January, 7 = July).",
        &Elem::U8(months),
        &lab,
    )?);
    let delta = leaps
        .iter()
        .map(|l| u8::try_from(l.2).map_err(|e| e.to_string()))
        .collect::<Result<_, _>>()?;
    out.push_str(&table(
        "_TM_LEAP_SECOND_DELTA_AT",
        "Column: TAI − UTC after the leap second, s.",
        &Elem::U8(delta),
        &lab,
    )?);

    // ── 1961–1971 segments ──
    out.push_str(&doc_comment(src, "const MJD_1972:")?);
    let _ = writeln!(out, "const _TM_MJD_1972 = {};\n", f64_lit(mjd_1972)?);
    let lab: Vec<String> = seg.iter().map(|s| format!("MJD {}", s[0])).collect();
    out.push_str(&doc_comment(src, "const TAI_UTC_PRE_1972:")?);
    let _ = writeln!(out, "const _TM_TAI_UTC_SEGMENT_COUNT = {};\n", seg.len());
    for (col, idx, doc) in [
        ("START", 0, "Column: MJD start."),
        ("OFFSET", 1, "Column: offset, s."),
        ("REFERENCE", 2, "Column: MJD reference."),
        ("RATE", 3, "Column: rate, s/day."),
    ] {
        let v = seg.iter().map(|s| s[idx]).collect();
        out.push_str(&table(
            &format!("_TM_TAI_UTC_SEGMENT_{col}"),
            doc,
            &Elem::F64(v),
            &lab,
        )?);
    }
    Ok(out.trim_end().to_string() + "\n")
}
