//! `libm.tsv` — reference outputs of every `libm` function tanmatra 1.3.0 calls
//! (`grep -oh 'libm::[a-z0-9_]*' rust-old/src`), for the M0 step 8 baselines.
//!
//! Each Cyrius counterpart (port plan §5) is run on these inputs, and the ulp distance
//! distribution per function and per target sets the class-C bound.
//!
//! Classes:
//! - **A**: correctly rounded or exact on both sides — sqrt, fabs, floor, trunc, round
//!   (the port's half-away helper, not `f64_round`), cbrt (CORE-MATH vs ganita), ldexp, and
//!   hypot (the port's `_tm_hypot` is Rust libm's algorithm operation for operation).
//! - **C**: everything else.
//!
//! Inputs are the special values, the ranges tanmatra feeds each function, and wide
//! seeded sweeps.

use crate::grids::{Rng, specials};
use crate::tsv::{Class, Module, f, i};

const PI: f64 = core::f64::consts::PI;

fn unary(m: &mut Module, name: &str, class: Class, fun: fn(f64) -> f64, xs: &[f64]) {
    for &x in xs {
        m.row(name, class, &[f(x)], &[f(fun(x))]);
    }
}

/// Specials, ±2^k at every 16th binade, and a wide signed log-uniform sweep.
fn wide(seed: u64, n: usize) -> Vec<f64> {
    let mut v = specials();
    for k in (-1074..=1023).step_by(16) {
        let x = libm::ldexp(1.0, k);
        v.push(x);
        v.push(-x);
    }
    let mut rng = Rng::new(seed);
    v.extend((0..n).map(|_| rng.signed_log_uniform(-300.0, 300.0)));
    v
}

fn sweep(seed: u64, n: usize, lo: f64, hi: f64) -> Vec<f64> {
    let mut rng = Rng::new(seed);
    (0..n).map(|_| rng.uniform(lo, hi)).collect()
}

fn log_sweep(seed: u64, n: usize, lo_exp10: f64, hi_exp10: f64) -> Vec<f64> {
    let mut rng = Rng::new(seed);
    (0..n)
        .map(|_| rng.log_uniform(lo_exp10, hi_exp10))
        .collect()
}

pub fn emit() -> Module {
    let mut m = Module::new("libm");

    // ── exact / correctly rounded (class A) ──
    let mut xs = wide(0x5157, 300);
    xs.extend(log_sweep(0x5158, 200, -10.0, 10.0));
    unary(&mut m, "libm::sqrt", Class::A, libm::sqrt, &xs);
    unary(
        &mut m,
        "libm::fabs",
        Class::A,
        libm::fabs,
        &wide(0xfab5, 100),
    );

    let mut halves: Vec<f64> = (-20..=20).map(|k| f64::from(k) + 0.5).collect();
    halves.extend((-20..=20).map(|k| f64::from(k) - 0.5 + f64::EPSILON));
    halves.extend([
        0.499_999_999_999_999_94,
        -0.499_999_999_999_999_94,
        4_503_599_627_370_495.5,
    ]);
    halves.extend(sweep(0x7a11, 200, -1e3, 1e3));
    halves.extend(wide(0x7a12, 100));
    unary(&mut m, "libm::floor", Class::A, libm::floor, &halves);
    unary(&mut m, "libm::trunc", Class::A, libm::trunc, &halves);
    unary(&mut m, "libm::round", Class::A, libm::round, &halves);

    // cbrt: every mass number, perfect cubes, both signs, and a wide sweep.
    let mut cb: Vec<f64> = (1..=300).map(f64::from).collect();
    cb.extend((1..=1000).map(|k| f64::from(k * k * k)));
    cb.extend((1..=100).map(|k| -f64::from(k * k * k)));
    cb.extend(wide(0xcb27, 400));
    unary(&mut m, "libm::cbrt", Class::A, libm::cbrt, &cb);

    let mut rng = Rng::new(0x1de8);
    let mut ld = vec![
        (1.0, 0),
        (1.0, -1074),
        (1.0, 1023),
        (1.5, -1075),
        (f64::MAX, 1),
        (f64::NAN, 3),
    ];
    for _ in 0..300 {
        #[allow(clippy::cast_possible_truncation)]
        let k = rng.int(-64, 64) as i32;
        ld.push((rng.signed_log_uniform(-20.0, 20.0), k));
    }
    for (x, k) in ld {
        m.row(
            "libm::ldexp",
            Class::A,
            &[f(x), i(k)],
            &[f(libm::ldexp(x, k))],
        );
    }

    // ── transcendental (class C) ──
    let mut tr = specials();
    tr.extend(sweep(0x5101, 400, -2.0 * PI, 2.0 * PI));
    tr.extend(sweep(0x5102, 200, 0.0, 1e3));
    tr.extend(log_sweep(0x5103, 100, -12.0, 15.0));
    unary(&mut m, "libm::sin", Class::C, libm::sin, &tr);
    unary(&mut m, "libm::cos", Class::C, libm::cos, &tr);
    unary(&mut m, "libm::tan", Class::C, libm::tan, &tr);

    let mut ex = specials();
    ex.extend(sweep(0xe001, 400, -50.0, 50.0));
    ex.extend(sweep(0xe002, 200, -745.2, 709.8));
    ex.extend(sweep(0xe003, 100, -1e-6, 1e-6));
    unary(&mut m, "libm::exp", Class::C, libm::exp, &ex);

    let mut lg = wide(0x1060, 300);
    lg.extend(sweep(0x1061, 200, 0.5, 2.0));
    lg.extend(sweep(0x1062, 100, 1.0 - 1e-9, 1.0 + 1e-9));
    unary(&mut m, "libm::log", Class::C, libm::log, &lg);

    let mut sh = specials();
    sh.extend(sweep(0x5a01, 300, -20.0, 20.0));
    sh.extend(sweep(0x5a02, 100, -710.0, 710.0));
    sh.extend(sweep(0x5a03, 100, -1e-4, 1e-4));
    unary(&mut m, "libm::sinh", Class::C, libm::sinh, &sh);

    // pow: tanmatra's exponents (2/3, 5/3, 1/3, half-integers, small integers) on its base
    // ranges, then random pairs and the special-value grid.
    let exps = [
        1.0 / 3.0,
        2.0 / 3.0,
        5.0 / 3.0,
        0.5,
        1.5,
        2.5,
        -0.5,
        -1.0,
        2.0,
        3.0,
        4.0,
        0.0,
    ];
    let bases = log_sweep(0x9001, 60, -3.0, 3.0);
    for &y in &exps {
        for &x in &bases {
            m.row("libm::pow", Class::C, &[f(x), f(y)], &[f(libm::pow(x, y))]);
        }
    }
    let mut rng = Rng::new(0x9002);
    for _ in 0..400 {
        let x = rng.log_uniform(-30.0, 30.0);
        let y = rng.uniform(-10.0, 10.0);
        m.row("libm::pow", Class::C, &[f(x), f(y)], &[f(libm::pow(x, y))]);
    }
    for &x in &specials() {
        for &y in &specials() {
            m.row("libm::pow", Class::C, &[f(x), f(y)], &[f(libm::pow(x, y))]);
        }
    }

    // hypot: tanmatra nests it for |p| (momenta 1e-9 .. 1e8 MeV), plus wide pairs.
    let mut rng = Rng::new(0x4970);
    for _ in 0..300 {
        let x = rng.signed_log_uniform(-9.0, 8.0);
        let y = rng.signed_log_uniform(-9.0, 8.0);
        m.row(
            "libm::hypot",
            Class::A,
            &[f(x), f(y)],
            &[f(libm::hypot(x, y))],
        );
    }
    for _ in 0..200 {
        let x = rng.signed_log_uniform(-300.0, 300.0);
        let y = rng.signed_log_uniform(-300.0, 300.0);
        m.row(
            "libm::hypot",
            Class::A,
            &[f(x), f(y)],
            &[f(libm::hypot(x, y))],
        );
    }
    for &x in &specials() {
        for &y in &specials() {
            m.row(
                "libm::hypot",
                Class::A,
                &[f(x), f(y)],
                &[f(libm::hypot(x, y))],
            );
        }
    }

    m
}
