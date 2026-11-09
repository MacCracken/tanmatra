//! `optics.tsv` — the port plan M6 grid for `tanmatra::optics` (the `optics` feature, against
//! the Rust prakash 1.2.0 the oracle links; the port runs against the Cyrius prakash 2.15.1).
//!
//! - **`hydrogen_balmer_series`** (class B): the six Storey & Hummer Case B lines.
//! - **`spectral_series_lines`** (class D: its intensities sum Einstein A-values): Z ∈ {0, 1, 2,
//!   26}, every series n_lower ≤ 5, upper limits from below the series start to n = 61.
//! - **`lines_to_spd`** (class C: exp per sample per line): the Balmer lines at FWHM 0.1, 1 and
//!   2 nm over 380–780 nm at 1 and 0.5 nm steps and over 650–662 nm at 0.1 nm; one line, no
//!   lines, Lyman lines outside the window; and every invalid-range guard (non-finite bounds,
//!   end < start, step ≤ 0 or non-finite, more than 10⁷ samples, FWHM ≤ 0 or NaN).
//! - **`line_to_rgb`** (class B): every branch edge of prakash's piecewise CIE fit, the Balmer
//!   lines, a 0.5 nm sweep over 380–780 nm, and the out-of-range values.

use crate::tsv::{Class, Module, Tok, f, i};
use tanmatra::optics::{self as op, SpectralLine};

fn lines_toks(lines: &[SpectralLine]) -> Vec<Tok> {
    let mut out = vec![Tok::V(lines.len())];
    for l in lines {
        out.extend([f(l.wavelength_nm), f(l.intensity)]);
    }
    out
}

fn spd_row(m: &mut Module, lines: &[SpectralLine], fwhm: f64, start: f64, end: f64, step: f64) {
    let spd = op::lines_to_spd(lines, fwhm, start, end, step);
    let mut ins = lines_toks(lines);
    ins.extend([f(fwhm), f(start), f(end), f(step)]);
    let mut out = vec![f(spd.start_nm), f(spd.step_nm), Tok::V(spd.values.len())];
    out.extend(spd.values.iter().map(|&v| f(v)));
    m.row("lines_to_spd", Class::C, &ins, &out);
}

pub fn emit() -> Module {
    let mut m = Module::new("optics");

    let balmer = op::hydrogen_balmer_series();
    m.row(
        "hydrogen_balmer_series",
        Class::B,
        &[],
        &lines_toks(&balmer),
    );

    for z in [0, 1, 2, 26] {
        for nl in 0..=5u32 {
            for nu in [0, nl, nl + 1, nl + 2, nl + 5, 12, 61] {
                m.row(
                    "spectral_series_lines",
                    Class::D,
                    &[i(z), i(nl), i(nu)],
                    &lines_toks(&op::spectral_series_lines(z, nl, nu)),
                );
            }
        }
    }

    // ── SPDs ──
    for fwhm in [0.1, 1.0, 2.0] {
        spd_row(&mut m, &balmer, fwhm, 380.0, 780.0, 1.0);
        spd_row(&mut m, &balmer, fwhm, 380.0, 780.0, 0.5);
        spd_row(&mut m, &balmer, fwhm, 650.0, 662.0, 0.1);
    }
    let one = [SpectralLine {
        wavelength_nm: 500.0,
        intensity: 2.5,
    }];
    let lyman = op::spectral_series_lines(1, 1, 6);
    spd_row(&mut m, &one, 1.0, 480.0, 520.0, 0.25);
    spd_row(&mut m, &[], 1.0, 380.0, 400.0, 1.0);
    spd_row(&mut m, &lyman, 1.0, 380.0, 400.0, 1.0);
    spd_row(&mut m, &lyman, 5.0, 90.0, 130.0, 0.5);
    // The guards.
    for (fwhm, start, end, step) in [
        (0.0, 380.0, 400.0, 1.0),
        (-1.0, 380.0, 400.0, 1.0),
        (f64::NAN, 380.0, 400.0, 1.0),
        (1.0, f64::NAN, 400.0, 1.0),
        (1.0, f64::INFINITY, 400.0, 1.0),
        (1.0, f64::NEG_INFINITY, 400.0, 1.0),
        (1.0, 380.0, f64::INFINITY, 1.0),
        (1.0, 380.0, f64::NAN, 1.0),
        (1.0, 400.0, 380.0, 1.0),
        (1.0, 380.0, 400.0, 0.0),
        (1.0, 380.0, 400.0, -1.0),
        (1.0, 380.0, 400.0, f64::NAN),
        (1.0, 380.0, 400.0, f64::INFINITY),
        (1.0, 0.0, 1e8, 1.0),
        (1.0, 0.0, 10_000_001.0, 1.0),
        (1.0, 0.0, 1e7, 1.0),
        (1.0, 380.0, 380.0, 1.0),
        (1.0, 380.0, 380.5, 1.0),
        (1.0, 380.0, 381.0, 0.3),
        (1.0, -0.0, 2.0, 1.0),
    ] {
        if end - start > 1e3 && step >= 1.0 && fwhm > 0.0 && (end - start) / step <= 1e7 {
            // 10⁷ samples is a valid SPD: too big for the corpus. Only the guard side is kept.
            continue;
        }
        spd_row(&mut m, &balmer[..2], fwhm, start, end, step);
    }

    // ── colours ──
    let mut wl: Vec<f64> = vec![
        379.999_999,
        380.0,
        400.0,
        419.999,
        420.0,
        439.999,
        440.0,
        465.0,
        489.999,
        490.0,
        500.0,
        509.999,
        510.0,
        545.0,
        579.999,
        580.0,
        600.0,
        644.999,
        645.0,
        656.3,
        700.0,
        700.001,
        740.0,
        780.0,
        780.000_001,
        0.0,
        -1.0,
        1e9,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ];
    wl.extend(balmer.iter().map(|l| l.wavelength_nm));
    wl.extend((0..=800).map(|k| 380.0 + 0.5 * f64::from(k)));
    for &w in &wl {
        let line = SpectralLine {
            wavelength_nm: w,
            intensity: 1.0,
        };
        let out = match op::line_to_rgb(&line) {
            Ok(c) => vec![f(c.r), f(c.g), f(c.b)],
            Err(e) => {
                let dbg = format!("{e:?}");
                let end = dbg.find(['(', ' ', '{']).unwrap_or(dbg.len());
                vec![Tok::Err(dbg[..end].to_string())]
            }
        };
        m.row("line_to_rgb", Class::B, &[f(w), f(1.0)], &out);
    }
    m
}
