//! `atomic.tsv` — the port plan M3 grid for `tanmatra::atomic`.
//!
//! - **Configurations and data (class A):**
//!   - every Z = 0..=120, and u32::MAX, through `electron_configuration`, its two formats, the
//!     ionization energy and the electron affinity;
//!   - hand-built configurations for the formats;
//!   - `QuantumNumbers::new` and the `OrbitalType` properties.
//! - **Hydrogen-like lines and levels (class B):**
//!   - every n1, n2 ≤ 12 (both orders, equal, zero) for Z ∈ {0, 1, 2, 3, 26, 92, 118, u32::MAX};
//!   - the reduced mass over the nucleus sample;
//!   - first-order fine structure for every (n, 2j) with n ≤ 12;
//!   - the five named series;
//!   - Landé g, the Zeeman and Stark shifts.
//! - **QED (class C where ln is used):** every (n, l, j) with n ≤ 4 for Z = 1..=118, plus the
//!   invalid levels; the vacuum polarization.
//! - **Radial functions (class C, exp):** every n ≤ 10, l < n at 50 radii, r = 0.1·k·n²/Z, for
//!   Z ∈ {1, 2, 92}; n = 20, 40 and 60; the non-finite radii.
//! - **Einstein A and B (class D, the radial dipole integral):** every E1-allowed
//!   (n_u ≤ 10, l_u) → (n_l, l_l) for Z ∈ {1, 2, 3, 26, 92}; a few n_u up to 60; the invalid
//!   transitions.
//! - **Dirac, hyperfine, g-factors, Breit (class B: sqrt and arithmetic only).**
//!
//! The oracle is built with overflow checks (scripts/oracle/Cargo.toml), so the grid keeps out
//! of the u32 sites that overflow in Rust: `2n − 1` with n ≥ 2³¹ (`hydrogen_level_energy_ev`,
//! `dirac_energy_mev`), `2l + 2` and `2j + 1` near u32::MAX (`lande_g_factor`,
//! `lamb_shift_nlj_ev`), `ml1 − ml2` beyond i32 (`check_selection_rules_full`), and series of
//! 4·10⁹ lines.

use crate::m_nucleus::tabulated;
use crate::tsv::{Class, Module, Tok, e, err, f, i, s};
use std::collections::BTreeSet;
use tanmatra::atomic::{self as at, ElectronAffinity, OrbitalFilling, OrbitalType, QuantumNumbers};
use tanmatra::error::TanmatraError;

/// Every OrbitalType variant, in declaration order.
pub const ORBITALS: [OrbitalType; 4] = [
    OrbitalType::S,
    OrbitalType::P,
    OrbitalType::D,
    OrbitalType::F,
];

const BIG: u32 = 1 << 30;

fn res_f(r: Result<f64, TanmatraError>) -> Vec<Tok> {
    match r {
        Ok(x) => vec![f(x)],
        Err(x) => vec![err(&x)],
    }
}

/// A configuration as tokens: `v:len`, then `i:n e:OrbitalType.X i:electrons` per subshell.
fn config_toks(c: &[OrbitalFilling]) -> Vec<Tok> {
    let mut out = vec![Tok::V(c.len())];
    for x in c {
        out.extend([i(x.n), e("OrbitalType", x.orbital), i(x.electrons)]);
    }
    out
}

fn fill(n: u32, orbital: OrbitalType, electrons: u32) -> OrbitalFilling {
    OrbitalFilling {
        n,
        orbital,
        electrons,
    }
}

fn series_toks(r: Result<Vec<(u32, f64)>, TanmatraError>) -> Vec<Tok> {
    match r {
        Ok(lines) => {
            let mut out = vec![Tok::V(lines.len())];
            for (n, nm) in lines {
                out.extend([i(n), f(nm)]);
            }
            out
        }
        Err(x) => vec![err(&x)],
    }
}

fn ea_toks(r: Result<ElectronAffinity, TanmatraError>) -> Vec<Tok> {
    match r {
        Ok(ElectronAffinity::Bound(ev)) => vec![Tok::E("ElectronAffinity", "Bound".into()), f(ev)],
        Ok(ElectronAffinity::Unbound) => vec![Tok::E("ElectronAffinity", "Unbound".into())],
        Ok(_) => vec![Tok::E("ElectronAffinity", "Unknown".into())],
        Err(x) => vec![err(&x)],
    }
}

/// Every (n, 2j) with n ≤ 4: the first-order fine-structure levels.
fn levels_n4() -> Vec<(u32, u32)> {
    let mut out = Vec::new();
    for n in 1..=4u32 {
        for tj in (1..2 * n).step_by(2) {
            out.push((n, tj));
        }
    }
    out
}

/// Every (n, l, 2j) with n ≤ 4, j = l ± 1/2.
fn nlj_n4() -> Vec<(u32, u32, u32)> {
    let mut out = Vec::new();
    for n in 1..=4u32 {
        for l in 0..n {
            if l > 0 {
                out.push((n, l, 2 * l - 1));
            }
            out.push((n, l, 2 * l + 1));
        }
    }
    out
}

fn configurations(m: &mut Module) {
    for o in ORBITALS {
        let v = [e("OrbitalType", o)];
        m.row("OrbitalType::l", Class::A, &v, &[i(o.l())]);
        m.row(
            "OrbitalType::max_electrons",
            Class::A,
            &v,
            &[i(o.max_electrons())],
        );
        m.row(
            "OrbitalType::symbol",
            Class::A,
            &v,
            &[s(o.symbol().to_string())],
        );
    }

    for n in [0, 1, 2, 3, 4, u32::MAX] {
        for l in [0, 1, 2, 3, 4, u32::MAX] {
            let mut mls: Vec<i32> = (-5..=5).collect();
            mls.extend([i32::MIN, i32::MAX]);
            for ml in mls {
                for ms in [-2, -1, 0, 1] {
                    let out = match QuantumNumbers::new(n, l, ml, ms) {
                        Ok(q) => vec![i(q.n), i(q.l), i(q.ml), i(q.ms)],
                        Err(x) => vec![err(&x)],
                    };
                    m.row(
                        "QuantumNumbers::new",
                        Class::A,
                        &[i(n), i(l), i(ml), i(ms)],
                        &out,
                    );
                }
            }
        }
    }

    let mut zs: Vec<u32> = (0..=120).collect();
    zs.push(u32::MAX);
    for &z in &zs {
        let cfg = at::electron_configuration(z);
        let out = match &cfg {
            Ok(c) => config_toks(c),
            Err(x) => vec![err(x)],
        };
        m.row("electron_configuration", Class::A, &[i(z)], &out);
        m.row(
            "ionization_energy_ev",
            Class::A,
            &[i(z)],
            &res_f(at::ionization_energy_ev(z)),
        );
        m.row(
            "electron_affinity",
            Class::A,
            &[i(z)],
            &ea_toks(at::electron_affinity(z)),
        );
        m.row(
            "electron_affinity_ev",
            Class::A,
            &[i(z)],
            &res_f(at::electron_affinity_ev(z)),
        );
    }

    // The formats, over every element's configuration and hand-built ones.
    let mut cfgs: Vec<Vec<OrbitalFilling>> = (1..=118)
        .filter_map(|z| at::electron_configuration(z).ok())
        .collect();
    cfgs.extend([
        vec![],
        vec![fill(0, OrbitalType::S, 0)],
        vec![fill(1, OrbitalType::S, 2), fill(2, OrbitalType::F, 14)],
        vec![fill(12, OrbitalType::D, 999), fill(3, OrbitalType::P, 0)],
        vec![fill(1, OrbitalType::S, 1); 25],
        vec![fill(BIG, OrbitalType::F, BIG)],
    ]);
    for c in &cfgs {
        let text = at::format_configuration(c);
        m.row(
            "format_configuration",
            Class::A,
            &config_toks(c),
            &[s(text)],
        );
    }
    // The noble-gas cores switch at Z > 2, 10, 18, 36, 54, 86; the configuration need not
    // be Z's own.
    let core_zs = [
        0, 1, 2, 3, 10, 11, 18, 19, 36, 37, 54, 55, 86, 87, 118, 200, BIG,
    ];
    for (k, c) in cfgs.iter().enumerate() {
        let own = u32::try_from(k + 1).unwrap_or(0);
        let mut zs: Vec<u32> = vec![own];
        if k % 9 == 0 || k >= 118 {
            zs.extend(core_zs);
        }
        for z in zs {
            let mut ins = config_toks(c);
            ins.push(i(z));
            let text = at::format_configuration_short(c, z);
            m.row("format_configuration_short", Class::A, &ins, &[s(text)]);
        }
    }
}

fn lines(m: &mut Module) {
    let zs = [0, 1, 2, 3, 26, 92, 118, u32::MAX];
    for z in zs {
        for n1 in 0..=12u32 {
            for n2 in 0..=12u32 {
                m.row(
                    "spectral_line_nm",
                    Class::B,
                    &[i(z), i(n1), i(n2)],
                    &res_f(at::spectral_line_nm(z, n1, n2)),
                );
            }
        }
        for (n1, n2) in [
            (1, u32::MAX),
            (u32::MAX, 1),
            (BIG, BIG + 1),
            (u32::MAX - 1, u32::MAX),
        ] {
            m.row(
                "spectral_line_nm",
                Class::B,
                &[i(z), i(n1), i(n2)],
                &res_f(at::spectral_line_nm(z, n1, n2)),
            );
        }
    }

    // The reduced mass: every tabulated nuclide, a sample of the rest, the edges.
    let mut nuclei: BTreeSet<(u32, u32)> = tabulated();
    nuclei.extend((1..=118u32).flat_map(|z| (z..=3 * z).step_by(11).map(move |a| (z, a))));
    nuclei.extend([
        (0, 0),
        (0, 1),
        (1, 0),
        (5, 4),
        (119, 300),
        (1, 1),
        (1, 2),
        (1, 3),
    ]);
    for &(z, a) in &nuclei {
        m.row(
            "reduced_mass_factor",
            Class::B,
            &[i(z), i(a)],
            &res_f(at::reduced_mass_factor(z, a)),
        );
    }
    let ions = [
        (1, 1),
        (1, 2),
        (1, 3),
        (2, 3),
        (2, 4),
        (3, 6),
        (3, 7),
        (26, 56),
        (92, 235),
        (92, 238),
        (118, 294),
        (0, 0),
        (1, 0),
    ];
    for (z, a) in ions {
        for n1 in 0..=6u32 {
            for n2 in 0..=6u32 {
                m.row(
                    "spectral_line_vacuum_nm",
                    Class::B,
                    &[i(z), i(a), i(n1), i(n2)],
                    &res_f(at::spectral_line_vacuum_nm(z, a, n1, n2)),
                );
            }
        }
    }

    // First-order fine structure.
    let zs = [0, 1, 2, 3, 26, 92, 118];
    for z in zs {
        for n in 0..=12u32 {
            for tj in 0..=2 * n + 1 {
                m.row(
                    "hydrogen_level_energy_ev",
                    Class::B,
                    &[i(z), i(n), i(tj)],
                    &res_f(at::hydrogen_level_energy_ev(z, n, tj)),
                );
            }
        }
        for tj in [1, 2 * BIG - 1, 2 * BIG + 1] {
            m.row(
                "hydrogen_level_energy_ev",
                Class::B,
                &[i(z), i(BIG), i(tj)],
                &res_f(at::hydrogen_level_energy_ev(z, BIG, tj)),
            );
        }
    }
    let mut lv = levels_n4();
    lv.extend([(0, 1), (1, 2), (1, 3)]);
    for z in [0, 1, 2, 3, 26, 92] {
        for &(n1, j1) in &lv {
            for &(n2, j2) in &lv {
                m.row(
                    "spectral_line_fine_nm",
                    Class::B,
                    &[i(z), i(n1), i(j1), i(n2), i(j2)],
                    &res_f(at::spectral_line_fine_nm(z, n1, j1, n2, j2)),
                );
            }
        }
    }

    type Series = fn(u32, u32) -> Result<Vec<(u32, f64)>, TanmatraError>;
    let series: [(&str, Series); 5] = [
        ("lyman_series", at::lyman_series),
        ("balmer_series", at::balmer_series),
        ("paschen_series", at::paschen_series),
        ("brackett_series", at::brackett_series),
        ("pfund_series", at::pfund_series),
    ];
    for (name, fun) in series {
        for z in [0, 1, 2, 26, 92] {
            for n_max in (0..=12).chain([20]) {
                m.row(
                    name,
                    Class::B,
                    &[i(z), i(n_max)],
                    &series_toks(fun(z, n_max)),
                );
            }
        }
    }
}

fn zeeman_stark(m: &mut Module) {
    let mut lj: Vec<(u32, u32)> = (0..=7)
        .flat_map(|l| (0..=16).map(move |tj| (l, tj)))
        .collect();
    lj.extend([(BIG, 2 * BIG + 1), (BIG, 2 * BIG - 1), (BIG, 1)]);
    for &(l, tj) in &lj {
        m.row(
            "lande_g_factor",
            Class::B,
            &[i(l), i(tj)],
            &[f(at::lande_g_factor(l, tj))],
        );
    }

    let fields = [
        0.0,
        1e-4,
        0.5,
        1.0,
        10.0,
        45.0,
        -2.0,
        f64::INFINITY,
        f64::NAN,
        f64::MAX,
    ];
    let mut zl: Vec<(u32, u32)> = Vec::new();
    for l in 0..=3u32 {
        if l > 0 {
            zl.push((l, 2 * l - 1));
        }
        zl.push((l, 2 * l + 1));
    }
    zl.extend([(0, 3), (2, 2)]);
    for &(l, tj) in &zl {
        for mj in [-7, -5, -3, -1, 0, 1, 3, 5, 7] {
            for &bt in &fields {
                m.row(
                    "zeeman_splitting_ev",
                    Class::B,
                    &[i(l), i(tj), i(mj), f(bt)],
                    &[f(at::zeeman_splitting_ev(l, tj, mj, bt))],
                );
            }
        }
    }
    for mj in [i32::MIN, i32::MAX] {
        m.row(
            "zeeman_splitting_ev",
            Class::B,
            &[i(1), i(3), i(mj), f(1.0)],
            &[f(at::zeeman_splitting_ev(1, 3, mj, 1.0))],
        );
    }

    let mut ks: Vec<i32> = (-6..=6).collect();
    ks.extend([i32::MIN, i32::MAX]);
    for n in (0..=6u32).chain([BIG]) {
        for &k in &ks {
            for ef in [0.0, 1.0, 5.142_206_747_63e11, -1e7, f64::INFINITY, f64::NAN] {
                m.row(
                    "stark_shift_hydrogen_ev",
                    Class::B,
                    &[i(n), i(k), f(ef)],
                    &[f(at::stark_shift_hydrogen_ev(n, k, ef))],
                );
            }
        }
    }
    for bt in [0.0, 1.0, 10.0, -3.0, 1e-5, f64::INFINITY, f64::NAN] {
        for ml in (-3..=3).chain([i32::MIN, i32::MAX]) {
            for ms in [-1, 1, 0, 3] {
                m.row(
                    "anomalous_zeeman_splitting_ev",
                    Class::B,
                    &[f(bt), i(ml), i(ms)],
                    &[f(at::anomalous_zeeman_splitting_ev(bt, ml, ms))],
                );
            }
        }
    }
}

fn qed(m: &mut Module) {
    let valid = nlj_n4();
    for z in 1..=118u32 {
        for &(n, l, tj) in &valid {
            m.row(
                "lamb_shift_nlj_ev",
                Class::C,
                &[i(z), i(n), i(l), i(tj)],
                &res_f(at::lamb_shift_nlj_ev(z, n, l, tj)),
            );
        }
    }
    for z in [0, 1, 92, u32::MAX] {
        for n in 0..=5u32 {
            for l in 0..=n {
                for tj in 0..=2 * n + 1 {
                    if z != 0 && z != u32::MAX && valid.contains(&(n, l, tj)) {
                        continue;
                    }
                    m.row(
                        "lamb_shift_nlj_ev",
                        Class::C,
                        &[i(z), i(n), i(l), i(tj)],
                        &res_f(at::lamb_shift_nlj_ev(z, n, l, tj)),
                    );
                }
            }
        }
        for (n, l, tj) in [(BIG, 1, 3), (3, BIG, 1), (4, 2, BIG)] {
            m.row(
                "lamb_shift_nlj_ev",
                Class::C,
                &[i(z), i(n), i(l), i(tj)],
                &res_f(at::lamb_shift_nlj_ev(z, n, l, tj)),
            );
        }
    }
    for z in (0..=118u32).chain([u32::MAX]) {
        for n in 0..=5u32 {
            for l in 0..=n {
                m.row(
                    "lamb_shift_ev",
                    Class::C,
                    &[i(z), i(n), i(l)],
                    &[f(at::lamb_shift_ev(z, n, l))],
                );
            }
        }
    }
    for z in [0, 1, 2, 3, 10, 26, 50, 82, 92, 118, u32::MAX] {
        for n in 0..=6u32 {
            for l in 0..=n {
                m.row(
                    "vacuum_polarization_ev",
                    Class::B,
                    &[i(z), i(n), i(l)],
                    &[f(at::vacuum_polarization_ev(z, n, l))],
                );
            }
        }
    }

    let dz = [
        0,
        1,
        2,
        3,
        10,
        26,
        50,
        82,
        92,
        100,
        118,
        136,
        137,
        138,
        140,
        1000,
        u32::MAX,
    ];
    for z in dz {
        for n in 0..=6u32 {
            for tj in 0..=2 * n + 1 {
                let ins = [i(z), i(n), i(tj)];
                m.row(
                    "dirac_energy_mev",
                    Class::B,
                    &ins,
                    &res_f(at::dirac_energy_mev(z, n, tj)),
                );
                m.row(
                    "dirac_binding_energy_ev",
                    Class::B,
                    &ins,
                    &res_f(at::dirac_binding_energy_ev(z, n, tj)),
                );
                m.row(
                    "relativistic_correction_ev",
                    Class::B,
                    &ins,
                    &res_f(at::relativistic_correction_ev(z, n, tj)),
                );
            }
        }
    }

    let gs = [
        5.585_694_689_3,
        0.857_438_233_5,
        -4.255_250_615,
        0.0,
        -0.0,
        f64::INFINITY,
        f64::NAN,
    ];
    for z in [0, 1, 2, 3, 26, 92] {
        for n in 0..=4u32 {
            for &g in &gs {
                m.row(
                    "hyperfine_splitting_ev",
                    Class::B,
                    &[i(z), i(n), f(g)],
                    &[f(at::hyperfine_splitting_ev(z, n, g))],
                );
            }
        }
    }
    let spins = [0.0, 0.5, 1.0, 1.5, 3.5, 4.5, -0.5, f64::INFINITY, f64::NAN];
    for z in [0, 1, 2, 26, 92] {
        for n in 0..=3u32 {
            for g in [5.585_694_689_3, 0.857_438_233_5, -0.0] {
                for &spin in &spins {
                    m.row(
                        "hyperfine_splitting_spin_ev",
                        Class::B,
                        &[i(z), i(n), f(g), f(spin)],
                        &[f(at::hyperfine_splitting_spin_ev(z, n, g, spin))],
                    );
                }
            }
        }
    }

    m.row(
        "electron_g_factor",
        Class::B,
        &[],
        &[f(at::electron_g_factor())],
    );
    for z in (0..=140u32).chain([1000, u32::MAX]) {
        m.row(
            "bound_electron_g_factor",
            Class::B,
            &[i(z)],
            &[f(at::bound_electron_g_factor(z))],
        );
        m.row(
            "breit_interaction_ev",
            Class::B,
            &[i(z)],
            &[f(at::breit_interaction_ev(z))],
        );
    }
}

fn radial(m: &mut Module) {
    // r is scaled by n²/Z, so ρ = 2Zr/n is nearly the same for every Z and only the
    // normalization differs: Z ∈ {1, 2, 92} covers it.
    let mut cases: Vec<(u32, u32, u32)> = Vec::new();
    for z in [1, 2, 92] {
        for n in 1..=10u32 {
            for l in 0..n {
                cases.push((z, n, l));
            }
        }
    }
    for z in [1, 92] {
        for n in [20, 40, 60u32] {
            for l in [0, n / 2, n - 1] {
                cases.push((z, n, l));
            }
        }
    }
    for &(z, n, l) in &cases {
        let scale = f64::from(n) * f64::from(n) / f64::from(z);
        for k in 0..50u32 {
            let r = 0.1 * f64::from(k) * scale;
            let ins = [i(z), i(n), i(l), f(r)];
            m.row(
                "radial_wavefunction",
                Class::C,
                &ins,
                &res_f(at::radial_wavefunction(z, n, l, r)),
            );
            m.row(
                "radial_probability_density",
                Class::C,
                &ins,
                &res_f(at::radial_probability_density(z, n, l, r)),
            );
        }
    }
    let odd_r = [
        -1.0,
        -0.0,
        f64::MIN_POSITIVE,
        1e300,
        f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ];
    let mut edges = vec![
        (1, 1, 0),
        (1, 3, 2),
        (92, 10, 3),
        (1, 60, 59),
        (u32::MAX, 2, 1),
    ];
    edges.extend([
        (1, 0, 0),
        (1, 61, 0),
        (1, 3, 3),
        (1, 3, BIG),
        (0, 1, 0),
        (0, 0, 0),
    ]);
    for &(z, n, l) in &edges {
        for r in odd_r.iter().copied().chain([1.0]) {
            let ins = [i(z), i(n), i(l), f(r)];
            m.row(
                "radial_wavefunction",
                Class::C,
                &ins,
                &res_f(at::radial_wavefunction(z, n, l, r)),
            );
            m.row(
                "radial_probability_density",
                Class::C,
                &ins,
                &res_f(at::radial_probability_density(z, n, l, r)),
            );
        }
    }
}

fn transitions(m: &mut Module) {
    let ls: Vec<u32> = (0..=6).chain([u32::MAX]).collect();
    for &l1 in &ls {
        for &l2 in &ls {
            let t = at::check_selection_rules(l1, l2);
            m.row(
                "check_selection_rules",
                Class::A,
                &[i(l1), i(l2)],
                &[e("TransitionType", t)],
            );
        }
    }
    // ±2^29: |ml1 − ml2| stays inside i32, where Rust's subtraction would overflow.
    let mls: Vec<i32> = (-3..=3)
        .chain([-(1 << 29), 1 << 29, (1 << 29) + 1])
        .collect();
    for l1 in 0..=3u32 {
        for &ml1 in &mls {
            for l2 in 0..=3u32 {
                for &ml2 in &mls {
                    let t = at::check_selection_rules_full(l1, ml1, l2, ml2);
                    m.row(
                        "check_selection_rules_full",
                        Class::A,
                        &[i(l1), i(ml1), i(l2), i(ml2)],
                        &[e("TransitionType", t)],
                    );
                }
            }
        }
    }

    // Every E1-allowed transition with n_u ≤ 10, a few long ones, and the invalid ones.
    let mut tr: Vec<(u32, u32, u32, u32)> = Vec::new();
    for nu in 2..=10u32 {
        for nl in 1..nu {
            for ll in 0..nl {
                if ll > 0 {
                    tr.push((nu, ll - 1, nl, ll));
                }
                tr.push((nu, ll + 1, nl, ll));
            }
        }
    }
    let mut zs_tr: Vec<(u32, (u32, u32, u32, u32))> = Vec::new();
    for z in [1, 2, 3, 26, 92] {
        zs_tr.extend(tr.iter().map(|&t| (z, t)));
    }
    let long = [
        (20, 1, 19, 0),
        (20, 19, 19, 18),
        (30, 29, 29, 28),
        (40, 20, 39, 19),
        (60, 1, 1, 0),
        (60, 59, 59, 58),
        (60, 30, 2, 1),
    ];
    zs_tr.extend(long.iter().map(|&t| (1, t)));
    let bad = [
        (3, 1, 2, 1),
        (3, 2, 2, 0),
        (2, 1, 3, 0),
        (2, 1, 2, 0),
        (2, 0, 0, 1),
        (3, 3, 2, 2),
        (3, 1, 2, 2),
        (61, 1, 1, 0),
        (2, 0, 1, 1),
    ];
    zs_tr.extend(bad.iter().map(|&t| (1, t)));
    zs_tr.extend([(0, (2, 1, 1, 0)), (u32::MAX, (2, 1, 1, 0))]);
    for (z, (nu, lu, nl, ll)) in zs_tr {
        let ins = [i(z), i(nu), i(lu), i(nl), i(ll)];
        m.row(
            "einstein_a_coefficient",
            Class::D,
            &ins,
            &res_f(at::einstein_a_coefficient(z, nu, lu, nl, ll)),
        );
        m.row(
            "einstein_b_coefficient",
            Class::D,
            &ins,
            &res_f(at::einstein_b_coefficient(z, nu, lu, nl, ll)),
        );
    }
}

pub fn emit() -> Module {
    let mut m = Module::new("atomic");
    configurations(&mut m);
    lines(&mut m);
    zeeman_stark(&mut m);
    qed(&mut m);
    radial(&mut m);
    transitions(&mut m);
    m
}
