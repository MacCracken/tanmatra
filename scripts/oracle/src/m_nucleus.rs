//! `nucleus.tsv` — the port plan M2 grid for `tanmatra::nucleus` (including the shell model).
//!
//! - **Every (Z, A) with 1 ≤ Z ≤ 118 and Z ≤ A ≤ 3Z** (14,160 nuclei) for `binding_energy`,
//!   `binding_energy_shell_corrected` and `ground_state_spin_parity`.
//! - **A sampled grid**, every 7th of those plus every tabulated nuclide and the presets, for the
//!   remaining `Nucleus` methods and the AME2020, radius and moment lookups.
//! - **Edges:** Nucleus::new's validation; nucleon counts 0..=220 for the magic-number, closure
//!   and occupation functions; every shell level, plus l > 6; every superallowed transition and
//!   non-table pairs.
//!
//! A nucleus input is the two tokens `i:Z i:A`. Classes: cbrt and sqrt are correctly rounded on
//! both sides, so the SEMF, the Myers–Swiatecki term, the masses and the radius are B (expected
//! 0 ulp). `atomic_mass_amu` goes through pow (the electron binding), so it is C. Tables, the
//! shell model and the predicates are A.

use crate::tsv::{Class, Module, Tok, b, err, f, i, s};
use std::collections::BTreeSet;
use tanmatra::nucleus::{
    Nucleus, ShellLevel, SuperallowedDecay, ground_state_spin_parity, is_magic_number,
    next_shell_closure, shell_closure_below, shell_model_levels, shell_occupation,
    superallowed_average_ft, superallowed_ft_values,
};

/// The full M2 nucleus grid: 1 ≤ Z ≤ 118, Z ≤ A ≤ 3Z.
pub fn grid() -> Vec<(u32, u32)> {
    (1..=118u32)
        .flat_map(|z| (z..=3 * z).map(move |a| (z, a)))
        .collect()
}

fn nuc(z: u32, a: u32) -> Nucleus {
    // Every caller passes a ≥ z ≥ 1, which `Nucleus::new` accepts.
    Nucleus::new(z, a).unwrap_or_else(|_| Nucleus::hydrogen_1())
}

fn zn(n: &Nucleus) -> [Tok; 2] {
    [i(n.z()), i(n.a())]
}

fn opt_f(x: Option<f64>) -> Vec<Tok> {
    x.map_or_else(|| vec![Tok::None], |v| vec![f(v)])
}

/// The tabulated nuclides (AME2020, radii, moments), from the oracle's lookups.
fn tabulated() -> BTreeSet<(u32, u32)> {
    let mut set = BTreeSet::new();
    for z in 1..=120 {
        for a in z..=300 {
            let n = nuc(z, a);
            if n.experimental_mass_excess_kev().is_some()
                || n.charge_radius_fm().is_some()
                || n.nuclear_moments().is_some()
            {
                set.insert((z, a));
            }
        }
    }
    set
}

/// A named `Nucleus` constructor.
type Preset = (&'static str, fn() -> Nucleus);

const PRESETS: [Preset; 6] = [
    ("Nucleus::hydrogen_1", Nucleus::hydrogen_1),
    ("Nucleus::helium_4", Nucleus::helium_4),
    ("Nucleus::carbon_12", Nucleus::carbon_12),
    ("Nucleus::iron_56", Nucleus::iron_56),
    ("Nucleus::uranium_235", Nucleus::uranium_235),
    ("Nucleus::uranium_238", Nucleus::uranium_238),
];

fn level_toks(l: &ShellLevel) -> [Tok; 3] {
    [i(l.n_shell), i(l.l), i(l.two_j)]
}

pub fn emit() -> Module {
    let mut m = Module::new("nucleus");

    // ── Nucleus::new validation ──
    let max = u32::MAX;
    for (z, a) in [
        (0, 0),
        (0, 1),
        (1, 0),
        (1, 1),
        (2, 1),
        (5, 4),
        (26, 56),
        (92, 238),
        (118, 354),
        (200, 500),
        (1, max),
        (max, max),
    ] {
        let out = match Nucleus::new(z, a) {
            Ok(n) => zn(&n).to_vec(),
            Err(e) => vec![err(&e)],
        };
        m.row("Nucleus::new", Class::A, &[i(z), i(a)], &out);
    }
    for (name, ctor) in PRESETS {
        m.row(name, Class::A, &[], &zn(&ctor()));
    }

    // ── the full grid ──
    let all = grid();
    for &(z, a) in &all {
        let n = nuc(z, a);
        let ins = zn(&n);
        m.row(
            "Nucleus::binding_energy",
            Class::B,
            &ins,
            &[f(n.binding_energy())],
        );
        let sc = n.binding_energy_shell_corrected();
        m.row(
            "Nucleus::binding_energy_shell_corrected",
            Class::B,
            &ins,
            &[f(sc)],
        );
        let (tj, par) = ground_state_spin_parity(&n);
        m.row("ground_state_spin_parity", Class::A, &ins, &[i(tj), i(par)]);
    }

    // ── sampled grid + every tabulated nuclide + the presets ──
    let mut sample: BTreeSet<(u32, u32)> = all.iter().copied().step_by(7).collect();
    sample.extend(tabulated());
    sample.extend(PRESETS.iter().map(|(_, c)| (c().z(), c().a())));
    for &(z, a) in &sample {
        let n = nuc(z, a);
        let ins = zn(&n);
        m.row("Nucleus::z", Class::A, &ins, &[i(n.z())]);
        m.row("Nucleus::a", Class::A, &ins, &[i(n.a())]);
        m.row("Nucleus::n", Class::A, &ins, &[i(n.n())]);
        let bpn = n.binding_energy_per_nucleon();
        m.row(
            "Nucleus::binding_energy_per_nucleon",
            Class::B,
            &ins,
            &[f(bpn)],
        );
        m.row(
            "Nucleus::mass_defect",
            Class::B,
            &ins,
            &[f(n.mass_defect())],
        );
        m.row(
            "Nucleus::nuclear_mass",
            Class::B,
            &ins,
            &[f(n.nuclear_mass())],
        );
        m.row(
            "Nucleus::atomic_mass_amu",
            Class::C,
            &ins,
            &[f(n.atomic_mass_amu())],
        );
        m.row(
            "Nucleus::nuclear_mass_amu",
            Class::B,
            &ins,
            &[f(n.nuclear_mass_amu())],
        );
        m.row(
            "Nucleus::nuclear_radius",
            Class::B,
            &ins,
            &[f(n.nuclear_radius())],
        );
        m.row("Nucleus::is_magic", Class::A, &ins, &[b(n.is_magic())]);
        m.row(
            "Nucleus::is_doubly_magic",
            Class::A,
            &ins,
            &[b(n.is_doubly_magic())],
        );
        let me = n.experimental_mass_excess_kev();
        m.row(
            "Nucleus::experimental_mass_excess_kev",
            Class::A,
            &ins,
            &opt_f(me),
        );
        let am = n.experimental_atomic_mass_amu();
        m.row(
            "Nucleus::experimental_atomic_mass_amu",
            Class::B,
            &ins,
            &opt_f(am),
        );
        m.row(
            "Nucleus::charge_radius_fm",
            Class::A,
            &ins,
            &opt_f(n.charge_radius_fm()),
        );
        let mom = n.nuclear_moments().map_or_else(
            || vec![Tok::None],
            |x| vec![f(x.magnetic_dipole_mu_n), f(x.electric_quadrupole_barn)],
        );
        m.row("Nucleus::nuclear_moments", Class::A, &ins, &mom);
    }

    // ── nucleon-count functions ──
    for count in 0..=220u32 {
        m.row(
            "is_magic_number",
            Class::A,
            &[i(count)],
            &[b(is_magic_number(count))],
        );
        m.row(
            "shell_closure_below",
            Class::A,
            &[i(count)],
            &[i(shell_closure_below(count))],
        );
        let nx = next_shell_closure(count).map_or(Tok::None, i);
        m.row("next_shell_closure", Class::A, &[i(count)], &[nx]);
        let occ = shell_occupation(count);
        let mut out = vec![Tok::V(occ.len())];
        for (lv, fill) in &occ {
            out.extend(level_toks(lv));
            out.push(i(*fill));
        }
        m.row("shell_occupation", Class::A, &[i(count)], &out);
    }

    // ── shell levels ──
    let levels = shell_model_levels();
    let mut out = vec![Tok::V(levels.len())];
    for &(n, l, j) in levels {
        out.extend([i(n), i(l), i(j)]);
    }
    m.row("shell_model_levels", Class::A, &[], &out);
    let mut probes: Vec<ShellLevel> = levels
        .iter()
        .map(|&(n_shell, l, two_j)| ShellLevel { n_shell, l, two_j })
        .collect();
    probes.extend(
        [(1, 7, 15), (2, 9, 17), (3, 0, 1)].map(|(n_shell, l, two_j)| ShellLevel {
            n_shell,
            l,
            two_j,
        }),
    );
    for lv in &probes {
        let ins = level_toks(lv);
        m.row(
            "ShellLevel::degeneracy",
            Class::A,
            &ins,
            &[i(lv.degeneracy())],
        );
        m.row("ShellLevel::j", Class::A, &ins, &[f(lv.j())]);
        m.row("ShellLevel::label", Class::A, &ins, &[s(lv.label())]);
    }

    // ── superallowed ──
    let sa = superallowed_ft_values();
    let mut out = vec![Tok::V(sa.len())];
    for d in &sa {
        out.extend(zn(&d.parent));
        out.extend(zn(&d.daughter));
        out.push(f(d.ft_seconds));
    }
    m.row("superallowed_ft_values", Class::A, &[], &out);
    let mut decays = sa.clone();
    decays.extend(
        [((8, 14), (7, 15)), ((26, 56), (25, 56)), ((7, 14), (8, 14))].map(
            |((pz, pa), (dz, da))| SuperallowedDecay {
                parent: nuc(pz, pa),
                daughter: nuc(dz, da),
                ft_seconds: 3000.0,
            },
        ),
    );
    for d in &decays {
        let mut ins = zn(&d.parent).to_vec();
        ins.extend(zn(&d.daughter));
        ins.push(f(d.ft_seconds));
        m.row(
            "SuperallowedDecay::corrected_ft_seconds",
            Class::A,
            &ins,
            &opt_f(d.corrected_ft_seconds()),
        );
    }
    let (avg, unc) = superallowed_average_ft();
    m.row("superallowed_average_ft", Class::A, &[], &[f(avg), f(unc)]);

    m
}
