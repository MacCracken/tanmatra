//! `reaction.tsv` — the port plan M4 grid for `tanmatra::reaction`.
//!
//! - **Data (class A):** every preset reaction and its two balance methods, plus hand-built
//!   reactions; the thermal cross sections, the resonance integrals over every tabulated nuclide
//!   and the nucleus sample, the fission yields, the two nucleosynthesis pathways and
//!   `CNO_NEUTRINO_LOSS_MEV`.
//! - **Arithmetic (class B):** `q_value` (with the signed-zero sums), the Coulomb barrier and the
//!   geometric cross section over nucleus pairs and the sample, Breit–Wigner around its
//!   resonance, `max_energy_loss_fraction` for A = 0..=300, and `moderating_ratio`.
//! - **ln:** `average_lethargy_gain` for A = 0..=300 and `collisions_to_thermalize`. They are class
//!   B: the port's `_tm_log` is Rust libm's log, bit for bit, as the port-parity rule requires here
//!   (ξ = 1 + (≈ −1) at large A puts K far past 64).

use crate::m_nucleus::tabulated;
use crate::tsv::{Class, Module, Tok, e, err, f, i, s};
use std::collections::BTreeSet;
use tanmatra::nucleus::Nucleus;
use tanmatra::reaction::{self as rx, NuclearReaction, NucleosynthesisPathway};

type Preset = (&'static str, fn() -> NuclearReaction);

/// Every preset reaction, in source order.
pub const PRESETS: [Preset; 7] = [
    ("dt_fusion", rx::dt_fusion),
    ("dd_fusion_he3", rx::dd_fusion_he3),
    ("dd_fusion_t", rx::dd_fusion_t),
    ("pp_chain_step1", rx::pp_chain_step1),
    ("u235_fission", rx::u235_fission),
    ("cno_cycle", rx::cno_cycle),
    ("triple_alpha", rx::triple_alpha),
];

type Pathway = (&'static str, fn() -> NucleosynthesisPathway);

/// Both nucleosynthesis pathways.
pub const PATHWAYS: [Pathway; 2] = [
    ("s_process_main", rx::s_process_main),
    ("r_process_main", rx::r_process_main),
];

fn nuc(z: u32, a: u32) -> Nucleus {
    Nucleus::new(z, a).unwrap_or_else(|_| Nucleus::hydrogen_1())
}

fn zn(n: &Nucleus) -> [Tok; 2] {
    [i(n.z()), i(n.a())]
}

fn opt_f(x: Option<f64>) -> Tok {
    x.map_or(Tok::None, f)
}

fn res_f(r: Result<f64, tanmatra::error::TanmatraError>) -> Vec<Tok> {
    match r {
        Ok(x) => vec![f(x)],
        Err(x) => vec![err(&x)],
    }
}

/// A reaction's fields as tokens: the name; the projectile (`n:` for None); the target; the
/// products; the neutron counts; the leptons; Q.
fn reaction_toks(r: &NuclearReaction) -> Vec<Tok> {
    let mut out = vec![s(r.name.clone())];
    match r.projectile {
        Some(p) => out.extend(zn(&p)),
        None => out.push(Tok::None),
    }
    out.extend(zn(&r.target));
    out.push(Tok::V(r.products.len()));
    for p in &r.products {
        out.extend(zn(p));
    }
    out.extend([i(r.neutrons_in), i(r.neutrons_out)]);
    out.push(Tok::V(r.leptons_out.len()));
    for l in &r.leptons_out {
        out.push(e("Lepton", *l));
    }
    out.push(f(r.q_value_mev));
    out
}

/// The inputs of a balance-method row: projectile (or `n:`), target, products, neutron counts.
fn balance_ins(r: &NuclearReaction) -> Vec<Tok> {
    let mut ins = Vec::new();
    match r.projectile {
        Some(p) => ins.extend(zn(&p)),
        None => ins.push(Tok::None),
    }
    ins.extend(zn(&r.target));
    ins.push(Tok::V(r.products.len()));
    for p in &r.products {
        ins.extend(zn(p));
    }
    ins.extend([i(r.neutrons_in), i(r.neutrons_out)]);
    ins
}

fn sums(m: &mut Module, rs: &[f64], ps: &[f64]) {
    let mut ins = vec![Tok::V(rs.len())];
    ins.extend(rs.iter().map(|&x| f(x)));
    ins.push(Tok::V(ps.len()));
    ins.extend(ps.iter().map(|&x| f(x)));
    m.row("q_value", Class::B, &ins, &[f(rx::q_value(rs, ps))]);
}

pub fn emit() -> Module {
    let mut m = Module::new("reaction");

    // ── presets and the balance methods ──
    let mut built: Vec<NuclearReaction> = Vec::new();
    for (name, ctor) in PRESETS {
        let r = ctor();
        m.row(name, Class::A, &[], &reaction_toks(&r));
        built.push(r);
    }
    let base = rx::dt_fusion();
    built.push(NuclearReaction {
        projectile: None,
        products: Vec::new(),
        neutrons_in: 7,
        ..base.clone()
    });
    built.push(NuclearReaction {
        projectile: Some(nuc(92, 238)),
        target: nuc(118, 294),
        products: vec![nuc(1, 1); 40],
        neutrons_out: u32::MAX,
        ..base.clone()
    });
    for r in &built {
        let ins = balance_ins(r);
        m.row(
            "NuclearReaction::nuclear_charge_change",
            Class::A,
            &ins,
            &[i(r.nuclear_charge_change())],
        );
        m.row(
            "NuclearReaction::baryon_number_change",
            Class::A,
            &ins,
            &[i(r.baryon_number_change())],
        );
    }
    m.row(
        "const",
        Class::A,
        &[s("CNO_NEUTRINO_LOSS_MEV")],
        &[f(rx::CNO_NEUTRINO_LOSS_MEV)],
    );

    // ── q_value: the f64 sums start from −0.0 in Rust 1.99 ──
    let me = tanmatra::constants::ELECTRON_MASS_MEV;
    let mp = tanmatra::constants::PROTON_MASS_MEV;
    let mn = tanmatra::constants::NEUTRON_MASS_MEV;
    for (rs, ps) in [
        (vec![], vec![]),
        (vec![-0.0], vec![]),
        (vec![], vec![-0.0]),
        (vec![-0.0], vec![0.0]),
        (vec![0.0], vec![-0.0]),
        (vec![-0.0, -0.0], vec![-0.0]),
        (vec![mp, mn], vec![1_875.612_945]),
        (
            vec![
                2.014_101_778 * 931.494_102_42,
                3.016_049_281 * 931.494_102_42,
            ],
            vec![4.002_603_254 * 931.494_102_42, mn],
        ),
        (vec![mp, mp], vec![1_875.612_945, me, 0.0]),
        (vec![1e300, 1e300], vec![1e300]),
        (vec![f64::MAX, f64::MAX], vec![]),
        (vec![f64::INFINITY], vec![f64::INFINITY]),
        (vec![f64::NAN, 1.0], vec![2.0]),
        (vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], vec![2.8]),
        (vec![1e16, 1.0, -1e16], vec![1.0]),
    ] {
        sums(&mut m, &rs, &ps);
    }

    // ── the nucleus sample: barrier pairs, geometric cross sections, resonance integrals ──
    let mut sample: BTreeSet<(u32, u32)> = tabulated();
    sample.extend((1..=118u32).flat_map(|z| (z..=3 * z).step_by(13).map(move |a| (z, a))));
    sample.extend([
        (1, 1),
        (1, 2),
        (1, 3),
        (2, 4),
        (92, 235),
        (92, 238),
        (94, 239),
        (5, 10),
    ]);
    for &(z, a) in &sample {
        let n = nuc(z, a);
        m.row(
            "geometric_cross_section_barns",
            Class::B,
            &zn(&n),
            &[f(rx::geometric_cross_section_barns(&n))],
        );
    }
    let pairs: Vec<(u32, u32)> = [
        (1, 1),
        (1, 2),
        (1, 3),
        (2, 3),
        (2, 4),
        (6, 12),
        (8, 16),
        (26, 56),
        (79, 197),
        (92, 238),
        (118, 294),
    ]
    .to_vec();
    for &(z1, a1) in &pairs {
        for &(z2, a2) in &pairs {
            let (n1, n2) = (nuc(z1, a1), nuc(z2, a2));
            let mut ins = zn(&n1).to_vec();
            ins.extend(zn(&n2));
            m.row(
                "coulomb_barrier",
                Class::B,
                &ins,
                &res_f(rx::coulomb_barrier(&n1, &n2)),
            );
        }
    }
    let mut ri: BTreeSet<(u32, u32)> = sample.clone();
    ri.extend([
        (0, 0),
        (0, 1),
        (40, 90),
        (42, 98),
        (50, 120),
        (53, 127),
        (55, 133),
        (u32::MAX, u32::MAX),
    ]);
    for &(z, a) in &ri {
        let ins = [i(z), i(a)];
        m.row(
            "resonance_integral_barns",
            Class::A,
            &ins,
            &[opt_f(rx::resonance_integral_barns(z, a))],
        );
        m.row(
            "fission_resonance_integral_barns",
            Class::A,
            &ins,
            &[opt_f(rx::fission_resonance_integral_barns(z, a))],
        );
        m.row(
            "n_alpha_resonance_integral_barns",
            Class::A,
            &ins,
            &[opt_f(rx::n_alpha_resonance_integral_barns(z, a))],
        );
    }

    // ── Breit–Wigner, around and away from the resonance ──
    let edges = [f64::NAN, f64::INFINITY, -0.0];
    for er in [0.1, 1.0] {
        for gam in [0.0, 1e-3, 0.1] {
            let mut es = vec![
                0.0,
                -1.0,
                1e-9,
                1e-6,
                er - gam,
                er - gam / 2.0,
                er,
                er + gam / 2.0,
                er + gam,
                10.0,
            ];
            es.extend(edges);
            for &en in &es {
                for (gi, gf) in [(gam, gam), (0.5 * gam, 0.25 * gam), (1e-6, 2.0)] {
                    for g in [0.5, 1.0] {
                        for mu in [0.0, 469.459, 938.272, -1.0] {
                            m.row(
                                "breit_wigner_cross_section",
                                Class::B,
                                &[f(en), f(er), f(gam), f(gi), f(gf), f(g), f(mu)],
                                &[f(rx::breit_wigner_cross_section(
                                    en, er, gam, gi, gf, g, mu,
                                ))],
                            );
                        }
                    }
                }
            }
        }
    }

    // ── tabulated data ──
    let th = rx::thermal_neutron_cross_sections();
    m.row(
        "thermal_neutron_cross_sections",
        Class::A,
        &[],
        &[Tok::V(th.len())],
    );
    for (k, x) in th.iter().enumerate() {
        let mut out = zn(&x.target).to_vec();
        out.extend([
            f(x.absorption_barns),
            f(x.fission_barns),
            f(x.scattering_barns),
        ]);
        m.row(
            "thermal_neutron_cross_sections[]",
            Class::A,
            &[i(u32::try_from(k).unwrap_or(0))],
            &out,
        );
    }
    for (name, ys) in [
        ("u235_fission_yields", rx::u235_fission_yields()),
        ("pu239_fission_yields", rx::pu239_fission_yields()),
    ] {
        m.row(name, Class::A, &[], &[Tok::V(ys.len())]);
        for (k, y) in ys.iter().enumerate() {
            m.row(
                &format!("{name}[]"),
                Class::A,
                &[i(u32::try_from(k).unwrap_or(0))],
                &[i(y.mass_number), f(y.yield_fraction)],
            );
        }
    }
    for (name, ctor) in PATHWAYS {
        let p = ctor();
        let mut out = vec![s(p.name.clone()), e("NucleosynthesisProcess", p.process)];
        out.extend(zn(&p.seed));
        out.push(Tok::V(p.steps.len()));
        for st in &p.steps {
            out.extend(zn(&st.nucleus));
            out.push(s(st.reaction.clone()));
        }
        m.row(name, Class::A, &[], &out);
    }

    // ── moderation ──
    for a in (0..=300u32).chain([1000, u32::MAX]) {
        m.row(
            "max_energy_loss_fraction",
            Class::B,
            &[i(a)],
            &[f(rx::max_energy_loss_fraction(a))],
        );
        m.row(
            "average_lethargy_gain",
            Class::B,
            &[i(a)],
            &[f(rx::average_lethargy_gain(a))],
        );
    }
    let eis = [2e6, 1e3, 1.0, 0.0253, 0.0, -1.0, f64::INFINITY, f64::NAN];
    let efs = [0.0253, 1.0, 2e6, 0.0, -1.0, f64::NAN];
    for a in [1, 2, 4, 9, 12, 16, 56, 238, 300] {
        for &ei in &eis {
            for &ef in &efs {
                m.row(
                    "collisions_to_thermalize",
                    Class::B,
                    &[i(a), f(ei), f(ef)],
                    &[f(rx::collisions_to_thermalize(a, ei, ef))],
                );
            }
        }
    }
    for xi in [1.0, 0.927, 0.158, 0.0, -0.5, f64::NAN] {
        for sc in [1.0, 3.45, 0.0, f64::INFINITY] {
            for ab in [0.022, 1e-5, 0.0, -1.0, f64::NAN] {
                m.row(
                    "moderating_ratio",
                    Class::B,
                    &[f(xi), f(sc), f(ab)],
                    &[f(rx::moderating_ratio(xi, sc, ab))],
                );
            }
        }
    }
    m
}
