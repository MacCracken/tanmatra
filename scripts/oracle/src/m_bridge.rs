//! `bridge.tsv` — the grid for `tanmatra::bridge`'s chrono bridge (port plan §6, the first
//! bridge release). The rest of `bridge.rs` waits on its consumers and has no rows yet.
//!
//! - **`utc_date_to_tai_offset`** (class A): years 1955–2035 × months 0–13, as M5's leap
//!   lookups, plus the i32 and u8 edges (`i32::MIN`, `i32::MAX`, year 0, month 255).
//! - **`tai_to_utc_seconds` and `utc_to_tai_seconds`** (class B: one add or subtract): every
//!   month where ΔAT changes and the month before it, found through the public
//!   `leap_seconds_at` and so covering every ΔAT value 0 and 10–37, plus the edge dates. Each is
//!   crossed with the float specials, the 2⁵³ neighbourhood (where ±ΔAT rounds), the Rust
//!   tests' 1e6 s, realistic TAI seconds since 1958 and a signed log-uniform sweep.

use crate::grids::{Rng, specials};
use crate::tsv::{Class, Module, f, i};
use tanmatra::bridge as br;
use tanmatra::timekeeping as tk;

/// (year, month) pairs: each month where ΔAT changes, the month before it, and the edges.
fn dates() -> Vec<(i32, u8)> {
    let mut v = Vec::new();
    let mut prev = tk::leap_seconds_at(1955, 12);
    for year in 1956..=2035 {
        for month in 1..=12u8 {
            let now = tk::leap_seconds_at(year, month);
            if now != prev {
                v.push(if month == 1 {
                    (year - 1, 12)
                } else {
                    (year, month - 1)
                });
                v.push((year, month));
            }
            prev = now;
        }
    }
    v.extend([
        (1960, 1),
        (2020, 1),
        (2020, 6),
        (2035, 12),
        (0, 0),
        (i32::MIN, 0),
        (i32::MAX, 255),
    ]);
    v
}

fn seconds() -> Vec<f64> {
    let mut v = specials();
    let two53 = 9_007_199_254_740_992.0;
    v.extend([
        0.5,
        0.1,
        37.0,
        -37.0,
        1e6,
        -1e6,
        1.9e9,
        2.1e9,
        1e-300,
        1e18,
        -1e18,
        two53,
        two53 - 1.0,
        two53 + 2.0,
        two53 / 2.0 + 0.5,
        -two53,
        f64::MIN,
    ]);
    let mut rng = Rng::new(0x6272_6964_6765);
    for _ in 0..8 {
        v.push(rng.uniform(0.0, 2.5e9));
    }
    for _ in 0..8 {
        v.push(rng.signed_log_uniform(-3.0, 20.0));
    }
    v
}

pub fn emit() -> Module {
    let mut m = Module::new("bridge");

    let mut ym: Vec<(i32, u8)> = Vec::new();
    for year in 1955..=2035 {
        for month in 0..=13u8 {
            ym.push((year, month));
        }
    }
    for year in [i32::MIN, -1, 0, 1, 9999, i32::MAX] {
        for month in [0, 1, 6, 12, 13, 255u8] {
            ym.push((year, month));
        }
    }
    for (year, month) in ym {
        m.row(
            "utc_date_to_tai_offset",
            Class::A,
            &[i(year), i(month)],
            &[i(br::utc_date_to_tai_offset(year, month))],
        );
    }

    let secs = seconds();
    for (year, month) in dates() {
        for &s in &secs {
            let ins = [f(s), i(year), i(month)];
            m.row(
                "tai_to_utc_seconds",
                Class::B,
                &ins,
                &[f(br::tai_to_utc_seconds(s, year, month))],
            );
            m.row(
                "utc_to_tai_seconds",
                Class::B,
                &ins,
                &[f(br::utc_to_tai_seconds(s, year, month))],
            );
        }
    }
    m
}
