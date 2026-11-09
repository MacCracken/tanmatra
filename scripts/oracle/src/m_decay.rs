//! `decay.tsv` — the port plan M2 grid for `tanmatra::decay`.
//!
//! - **Every known isotope:** its fields (`known_isotopes[]`, one row per element: input `i:k`);
//!   its half-life through `decay_constant`, `activity_bq`, and `remaining_fraction` at
//!   multiples of it.
//! - **Invalid-input edges:** 0, −0, negative, subnormal, ±∞ and NaN half-lives.
//! - **α, β⁻ and β⁺ decay** on the nucleus sample grid plus the edges.
//! - **Decay chains** from every isotope's nucleus, U-238, U-235, Th-232 and non-isotopes, at
//!   several step limits.
//! - **Bateman** (class D):
//!   - empty and single-member chains;
//!   - non-positive and non-finite times; negative, zero and non-finite λ;
//!   - equal and nearly equal λ; a 9-member long-lived chain;
//!   - the U-238 series at 1 y, 10³ y and 4.463e9 y; the U-235 and Th-232 series;
//!   - λ spread over 21 decades.
//!
//!   The series λ are ln 2 / T½ of the oracle's own `decay_chain` and isotope table.
//!
//! Classes: ln 2 / T and λN are B; `remaining_fraction` is pow, so C; Bateman is D.

use crate::m_nucleus::grid;
use crate::tsv::{Class, Module, Tok, b, e, err, f, i, s};
use std::collections::BTreeSet;
use tanmatra::decay::{
    DecayMode, activity_bq, alpha_decay, bateman_chain, beta_minus_decay, beta_plus_decay,
    decay_chain, decay_constant, known_isotopes, remaining_fraction,
};
use tanmatra::nucleus::Nucleus;

/// Every DecayMode variant, in declaration order.
pub const MODES: [DecayMode; 9] = [
    DecayMode::Alpha,
    DecayMode::BetaMinus,
    DecayMode::BetaPlus,
    DecayMode::ElectronCapture,
    DecayMode::Gamma,
    DecayMode::SpontaneousFission,
    DecayMode::ProtonEmission,
    DecayMode::NeutronEmission,
    DecayMode::IsomericTransition,
];

fn nuc(z: u32, a: u32) -> Nucleus {
    Nucleus::new(z, a).unwrap_or_else(|_| Nucleus::hydrogen_1())
}

fn zn(n: &Nucleus) -> [Tok; 2] {
    [i(n.z()), i(n.a())]
}

fn res_f(r: Result<f64, tanmatra::error::TanmatraError>) -> Vec<Tok> {
    match r {
        Ok(x) => vec![f(x)],
        Err(x) => vec![err(&x)],
    }
}

fn res_n(r: Result<Nucleus, tanmatra::error::TanmatraError>) -> Vec<Tok> {
    match r {
        Ok(n) => zn(&n).to_vec(),
        Err(x) => vec![err(&x)],
    }
}

/// The half-lives (s) of the nuclides `decay_chain` walks from `start`, by the same isotope
/// rule it uses, with the stable end member as λ = 0.
fn series_lambdas(start: &Nucleus) -> Vec<f64> {
    let all = known_isotopes();
    let chain = decay_chain(start, 50);
    let mut lam = Vec::new();
    let mut in_ground = false;
    for (n, mode) in &chain {
        let iso = all
            .iter()
            .find(|x| x.nucleus == *n && !(in_ground && x.is_isomer));
        if let Some(x) = iso {
            lam.push(decay_constant(x.half_life_seconds).unwrap_or(0.0));
        }
        in_ground = matches!(mode, DecayMode::IsomericTransition | DecayMode::Gamma);
    }
    lam.push(0.0);
    lam
}

fn bateman_row(m: &mut Module, lam: &[f64], n0: f64, t: f64) {
    let mut ins = vec![Tok::V(lam.len())];
    ins.extend(lam.iter().map(|&x| f(x)));
    ins.extend([f(n0), f(t)]);
    let pops = bateman_chain(lam, n0, t);
    let mut out = vec![Tok::V(pops.len())];
    out.extend(pops.iter().map(|&x| f(x)));
    m.row("bateman_chain", Class::D, &ins, &out);
}

pub fn emit() -> Module {
    let mut m = Module::new("decay");
    let isos = known_isotopes();

    // ── the isotope table ──
    m.row("known_isotopes", Class::A, &[], &[Tok::V(isos.len())]);
    for (k, x) in isos.iter().enumerate() {
        let mut out = zn(&x.nucleus).to_vec();
        out.extend([
            s(x.name.clone()),
            f(x.half_life_seconds),
            e("DecayMode", x.primary_decay),
            b(x.is_isomer),
            f(x.excitation_energy_kev),
        ]);
        m.row(
            "known_isotopes[]",
            Class::A,
            &[i(u32::try_from(k).unwrap_or(0))],
            &out,
        );
    }

    // ── half-life functions ──
    let mut hls: Vec<f64> = isos.iter().map(|x| x.half_life_seconds).collect();
    hls.extend([
        0.0,
        -0.0,
        -1.0,
        f64::from_bits(1),
        1e-300,
        1.0,
        f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ]);
    for &hl in &hls {
        m.row(
            "decay_constant",
            Class::B,
            &[f(hl)],
            &res_f(decay_constant(hl)),
        );
        for n_atoms in [0.0, 1.0, 6.022_140_76e23, -1.0] {
            m.row(
                "activity_bq",
                Class::B,
                &[f(hl), f(n_atoms)],
                &res_f(activity_bq(hl, n_atoms)),
            );
        }
        let times = if hl.is_finite() && hl > 0.0 {
            vec![
                0.0,
                0.5 * hl,
                hl,
                3.7 * hl,
                100.0 * hl,
                1e4 * hl,
                -1.0,
                f64::INFINITY,
                f64::NAN,
            ]
        } else {
            vec![0.0, 1.0, -1.0]
        };
        for t in times {
            m.row(
                "remaining_fraction",
                Class::C,
                &[f(hl), f(t)],
                &[f(remaining_fraction(hl, t))],
            );
        }
    }

    // ── α, β⁻, β⁺ ──
    let mut sample: BTreeSet<(u32, u32)> = grid().into_iter().step_by(7).collect();
    sample.extend([
        (1, 1),
        (1, 2),
        (2, 3),
        (2, 4),
        (3, 3),
        (3, 4),
        (3, 5),
        (5, 5),
        (5, 6),
        (92, 238),
        (118, 354),
    ]);
    for &(z, a) in &sample {
        let n = nuc(z, a);
        let ins = zn(&n);
        m.row("alpha_decay", Class::A, &ins, &res_n(alpha_decay(&n)));
        m.row(
            "beta_minus_decay",
            Class::A,
            &ins,
            &res_n(beta_minus_decay(&n)),
        );
        m.row(
            "beta_plus_decay",
            Class::A,
            &ins,
            &res_n(beta_plus_decay(&n)),
        );
    }

    // ── decay chains ──
    let mut starts: Vec<(u32, u32)> = Vec::new();
    for x in &isos {
        let za = (x.nucleus.z(), x.nucleus.a());
        if !starts.contains(&za) {
            starts.push(za);
        }
    }
    starts.extend([(26, 56), (1, 1), (82, 208)]);
    for &(z, a) in &starts {
        let n = nuc(z, a);
        for max_steps in [0usize, 1, 3, 10, 50] {
            let chain = decay_chain(&n, max_steps);
            let mut ins = zn(&n).to_vec();
            ins.push(i(u32::try_from(max_steps).unwrap_or(0)));
            let mut out = vec![Tok::V(chain.len())];
            for (d, mode) in &chain {
                out.extend(zn(d));
                out.push(e("DecayMode", *mode));
            }
            m.row("decay_chain", Class::A, &ins, &out);
        }
    }

    // ── Bateman ──
    let y = 365.2422 * 86400.0;
    bateman_row(&mut m, &[], 1.0, 1.0);
    for t in [0.0, -1.0, f64::NAN, f64::INFINITY, 1.0, 1e10] {
        bateman_row(&mut m, &[1e-3], 1e6, t);
    }
    bateman_row(&mut m, &[0.1, 0.05, 0.0], 1e6, 100.0);
    for t in [1.0, 1e3, 1e5] {
        bateman_row(&mut m, &[1e-3, 1e-3, 1e-3, 0.0], 1.0, t);
        bateman_row(&mut m, &[1e-3, 1e-3 * (1.0 + 1e-10), 2e-3, 0.0], 1.0, t);
    }
    bateman_row(
        &mut m,
        &[-1.0, 0.1, f64::INFINITY, f64::NAN, 0.0],
        1.0,
        10.0,
    );
    bateman_row(&mut m, &[1e-18, 1e3, 1e-5, 0.0], 6.022_140_76e23, 1e10);
    let u238 = series_lambdas(&Nucleus::uranium_238());
    for t in [y, 1e3 * y, 4.463e9 * y] {
        bateman_row(&mut m, &u238, 1.0, t);
    }
    // The 9-member long-lived chain: the first nine U-238 series members.
    let nine: Vec<f64> = u238.iter().copied().take(9).collect();
    for t in [1e3 * y, 1e6 * y, 1e9 * y] {
        bateman_row(&mut m, &nine, 1.0, t);
    }
    for start in [Nucleus::uranium_235(), nuc(90, 232)] {
        let lam = series_lambdas(&start);
        for t in [y, 1e6 * y, 1e10 * y] {
            bateman_row(&mut m, &lam, 1.0, t);
        }
    }

    m
}
