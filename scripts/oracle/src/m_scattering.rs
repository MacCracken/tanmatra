//! `scattering.tsv` — the port plan M4 grid for `tanmatra::scattering`.
//!
//! - **Angles:** θ ∈ {0, 1e-8, 1e-3, 0.1, 0.5, 1, π/2, 2, 3, π, −0.5, 2π, ±∞, NaN}, a superset of
//!   the plan's {1e-3, 0.1, π/2, π}.
//! - **Branch points, both sides:** Klein–Nishina's series below γ = 1e-2 (and the doubles
//!   adjacent to it), the form factor's series below |qR| = 0.5, pair production's k = 2
//!   threshold and its k = 4 join.
//! - **Classes:** the port reproduces Rust libm's sin, cos and log bit for bit (`_tm_sin`,
//!   `_tm_cos`, `_tm_log`), as the port-parity rule requires for formulas this ill-conditioned
//!   (K > 64 by perturbation: the Mott factor at β → 1, the form factor near its zeros,
//!   Klein–Nishina's closed form above γ = 1e-2). So every function built from them, sqrt, cbrt
//!   and arithmetic is class B. Only `rutherford_total_above_angle` (ganita's tan) is C.

use crate::grids::log_grid;
use crate::tsv::{Class, Module, Tok, f, i};
use tanmatra::constants::ELECTRON_MASS_MEV;
use tanmatra::scattering as sc;

use core::f64::consts::{FRAC_PI_2, PI};

fn angles() -> Vec<f64> {
    vec![
        0.0,
        1e-8,
        1e-3,
        0.1,
        0.5,
        1.0,
        FRAC_PI_2,
        2.0,
        3.0,
        PI,
        -0.5,
        2.0 * PI,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ]
}

/// x and the doubles on either side of it.
fn around(x: f64) -> [f64; 3] {
    [
        f64::from_bits(x.to_bits() - 1),
        x,
        f64::from_bits(x.to_bits() + 1),
    ]
}

fn shifts_toks(d: &[f64]) -> Vec<Tok> {
    let mut out = vec![Tok::V(d.len())];
    out.extend(d.iter().map(|&x| f(x)));
    out
}

pub fn emit() -> Module {
    let mut m = Module::new("scattering");
    let th = angles();
    let energies = [
        0.0,
        -1.0,
        1e-3,
        0.1,
        1.0,
        5.0,
        10.0,
        100.0,
        1e4,
        f64::INFINITY,
        f64::NAN,
    ];
    let zpairs = [(1u32, 1u32), (2, 79), (6, 6), (0, 79), (92, 92), (1, 82)];

    // ── Rutherford, Mott factor, closest approach, Sommerfeld ──
    for &(z1, z2) in &zpairs {
        for &en in &energies {
            for &t in &th {
                let ins = [i(z1), i(z2), f(en), f(t)];
                m.row(
                    "rutherford_differential",
                    Class::B,
                    &ins,
                    &[f(sc::rutherford_differential(z1, z2, en, t))],
                );
                m.row(
                    "rutherford_total_above_angle",
                    Class::C,
                    &ins,
                    &[f(sc::rutherford_total_above_angle(z1, z2, en, t))],
                );
            }
            m.row(
                "distance_of_closest_approach",
                Class::B,
                &[i(z1), i(z2), f(en)],
                &[f(sc::distance_of_closest_approach(z1, z2, en))],
            );
        }
        for beta in [0.0, 0.05, 0.5, 0.999, 1.0, -0.1, f64::NAN] {
            m.row(
                "sommerfeld_parameter",
                Class::B,
                &[i(z1), i(z2), f(beta)],
                &[f(sc::sommerfeld_parameter(z1, z2, beta))],
            );
        }
    }
    for beta in [0.0, 0.01, 0.5, 0.9, 0.999, 1.0, 1.5, -0.5, f64::NAN] {
        for &t in &th {
            m.row(
                "mott_correction_factor",
                Class::B,
                &[f(beta), f(t)],
                &[f(sc::mott_correction_factor(beta, t))],
            );
        }
    }

    // ── Legendre and partial waves ──
    let xs = [
        -1.0,
        -0.7,
        -0.5,
        -0.1,
        0.0,
        0.3,
        0.5,
        0.6,
        0.7,
        0.99,
        1.0,
        1.5,
        -2.0,
        f64::NAN,
    ];
    for l in (0..=20u32).chain([50]) {
        for &x in &xs {
            m.row(
                "legendre_polynomial",
                Class::B,
                &[i(l), f(x)],
                &[f(sc::legendre_polynomial(l, x))],
            );
        }
    }
    let sets: Vec<Vec<f64>> = vec![
        vec![],
        vec![0.5],
        vec![FRAC_PI_2],
        vec![0.5, 0.2, -0.1],
        vec![1e-8, 1e-8],
        vec![FRAC_PI_2; 6],
        (0..11).map(|l| 0.8 / f64::from(l + 1)).collect(),
        vec![3.0, -2.0, PI],
    ];
    for k in [0.0, 0.1, 1.0, 2.5, -1.0, f64::NAN] {
        for d in &sets {
            let mut ins = vec![f(k)];
            ins.extend(shifts_toks(d));
            m.row(
                "partial_wave_cross_section",
                Class::B,
                &ins,
                &[f(sc::partial_wave_cross_section(k, d))],
            );
            for &t in &[0.0, 0.1, 1.0, FRAC_PI_2, 2.0, PI, f64::NAN] {
                let mut ins2 = ins.clone();
                ins2.push(f(t));
                m.row(
                    "partial_wave_differential",
                    Class::B,
                    &ins2,
                    &[f(sc::partial_wave_differential(k, d, t))],
                );
            }
        }
    }

    // ── Born, screened Coulomb ──
    for z in (0..=120u32).chain([u32::MAX]) {
        m.row(
            "thomas_fermi_screening_fm",
            Class::B,
            &[i(z)],
            &[f(sc::thomas_fermi_screening_fm(z))],
        );
    }
    let born_th = [0.0, 1e-3, 0.1, 1.0, FRAC_PI_2, PI, f64::NAN];
    for &(z1, z2) in &[(1u32, 1u32), (2, 79), (6, 6), (1, 82)] {
        for en in [0.0, 0.01, 1.0, 5.0, 100.0, f64::NAN] {
            for &t in &born_th {
                for a in [0.0, 100.0, 30_000.0, 1e9, -1.0] {
                    m.row(
                        "born_screened_coulomb",
                        Class::B,
                        &[i(z1), i(z2), f(en), f(t), f(a)],
                        &[f(sc::born_screened_coulomb(z1, z2, en, t, a))],
                    );
                    for (m1, m2) in [
                        (938.272, 183_473.0),
                        (3727.379, 3727.379),
                        (0.511, 1e9),
                        (0.0, 1.0),
                    ] {
                        m.row(
                            "born_screened_coulomb_with_masses",
                            Class::B,
                            &[i(z1), i(z2), f(m1), f(m2), f(en), f(t), f(a)],
                            &[f(sc::born_screened_coulomb_with_masses(
                                z1, z2, m1, m2, en, t, a,
                            ))],
                        );
                    }
                }
            }
        }
    }

    // ── Mott for electrons, the form factor ──
    let te = [0.0, -1.0, 1e-3, 0.1, 1.0, 10.0, 100.0, 500.0, 1e4, f64::NAN];
    for z in [1u32, 6, 20, 79, 82] {
        for &en in &te {
            for &t in &th {
                m.row(
                    "mott_electron_differential",
                    Class::B,
                    &[i(z), f(en), f(t)],
                    &[f(sc::mott_electron_differential(z, en, t))],
                );
            }
        }
    }
    for (z, a) in [(1u32, 1u32), (6, 12), (20, 40), (79, 197), (82, 208)] {
        for &en in &te {
            for &t in &th {
                m.row(
                    "mott_electron_with_form_factor",
                    Class::B,
                    &[i(z), i(a), f(en), f(t)],
                    &[f(sc::mott_electron_with_form_factor(z, a, en, t))],
                );
            }
        }
    }
    // qR on both sides of 0.5 (the series switch), the first zeros, large and odd values.
    let mut qr: Vec<f64> = vec![
        0.0, -0.0, 1e-8, 0.1, 0.3, 0.49, 0.51, 0.7, 1.0, 2.0, 4.0, 4.4934, 7.7253, 10.0, 25.0,
        100.0, -0.3, -3.0,
    ];
    qr.extend(around(0.5));
    qr.extend(around(-0.5));
    qr.extend([f64::INFINITY, f64::NAN]);
    for &x in &qr {
        for r in [1.0, 2.5, 7.0] {
            let q = x / r;
            m.row(
                "nuclear_form_factor_uniform",
                Class::B,
                &[f(q), f(r)],
                &[f(sc::nuclear_form_factor_uniform(q, r))],
            );
        }
    }

    // ── Compton: the ratio, Klein–Nishina both sides of γ = 1e-2 ──
    let mut eg = log_grid(1e-6, 10, 4);
    eg.extend(around(1e-2 * ELECTRON_MASS_MEV));
    eg.extend([0.0, -1.0, f64::INFINITY, f64::NAN]);
    for &en in &eg {
        m.row(
            "klein_nishina_total",
            Class::B,
            &[f(en)],
            &[f(sc::klein_nishina_total(en))],
        );
        for &t in &th {
            m.row(
                "compton_energy_ratio",
                Class::B,
                &[f(en), f(t)],
                &[f(sc::compton_energy_ratio(en, t))],
            );
            m.row(
                "klein_nishina_differential",
                Class::B,
                &[f(en), f(t)],
                &[f(sc::klein_nishina_differential(en, t))],
            );
        }
    }
    // γ just below and above 1e-2 in γ itself (the branch reads γ = E/m_e).
    for g in [0.009_999_999, 0.01, 0.010_000_001, 0.005, 0.02] {
        let en = g * ELECTRON_MASS_MEV;
        m.row(
            "klein_nishina_total",
            Class::B,
            &[f(en)],
            &[f(sc::klein_nishina_total(en))],
        );
    }

    // ── pair production: k on both sides of 2 and 4 ──
    let mut ek = log_grid(1.0, 6, 8);
    for k in [2.0, 3.0, 4.0] {
        for x in around(k * ELECTRON_MASS_MEV) {
            ek.push(x);
        }
    }
    ek.extend([
        2.5 * ELECTRON_MASS_MEV,
        3.999 * ELECTRON_MASS_MEV,
        4.001 * ELECTRON_MASS_MEV,
    ]);
    ek.extend([0.0, -1.0, f64::INFINITY, f64::NAN, 1e300]);
    for z in [0u32, 1, 6, 26, 82, 92] {
        for &en in &ek {
            m.row(
                "pair_production_cross_section",
                Class::B,
                &[i(z), f(en)],
                &[f(sc::pair_production_cross_section(z, en))],
            );
        }
    }
    m
}
