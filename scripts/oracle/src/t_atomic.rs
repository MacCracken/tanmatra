//! `src/atomic_data.cyr` — the atomic tables (port plan M3).
//!
//! - **The enums and `OrbitalType`'s properties:** through the oracle's API.
//! - **Ionization energies and electron affinities (private tables):** through the public
//!   `ionization_energy_ev` and `electron_affinity` for Z = 1..=118. Each entry's label is its
//!   citation, the comment on its line in rust-old/src/atomic.rs.
//! - **The Aufbau filling order (private):** extracted by `rust_src`.
//! - **The 20 NIST ground-state exceptions (private logic):** DERIVED from the public
//!   `electron_configuration`, as the difference from the pure Madelung filling. The order and
//!   the exceptions are proven together by recomputing `electron_configuration` for
//!   Z = 0..=120 and u32::MAX.
//! - **The Bethe logarithms (private):** extracted from `bethe_log`'s match, then proven by
//!   recomputing `lamb_shift_nlj_ev` for Z = 0..=118 and u32::MAX over every (n ≤ 6, l ≤ n, 2j).

use crate::cyr::{Elem, doc_comment, header, snake, table};
use crate::m_atomic::ORBITALS;
use crate::rust_src::{Arm, match_f64_arms, u32_variant_pairs, z_comments};
use std::fmt::Write as _;
use tanmatra::atomic::{self as at, ElectronAffinity, OrbitalFilling, OrbitalType};
use tanmatra::constants::{ELECTRON_MASS_MEV, FINE_STRUCTURE};

fn orbital(name: &str) -> Result<OrbitalType, String> {
    ORBITALS
        .iter()
        .copied()
        .find(|o| format!("{o:?}") == name)
        .ok_or_else(|| format!("unknown OrbitalType::{name}"))
}

/// The pure Madelung filling of `z` electrons: the oracle's loop, without the exceptions.
fn madelung(order: &[(u32, OrbitalType)], z: u32) -> Vec<OrbitalFilling> {
    let mut config = Vec::new();
    let mut remaining = z;
    for &(n, orbital) in order {
        if remaining == 0 {
            break;
        }
        let fill = remaining.min(orbital.max_electrons());
        config.push(OrbitalFilling {
            n,
            orbital,
            electrons: fill,
        });
        remaining -= fill;
    }
    config
}

/// One exception step: set subshell (n, orbital) to `electrons`.
type Op = (u32, OrbitalType, u32);

/// The oracle's `set_subshell`, transcribed.
fn set_subshell(config: &mut Vec<OrbitalFilling>, (n, orbital, electrons): Op) {
    for entry in config.iter_mut() {
        if entry.n == n && entry.orbital == orbital {
            entry.electrons = electrons;
            return;
        }
    }
    if electrons > 0 {
        config.push(OrbitalFilling {
            n,
            orbital,
            electrons,
        });
    }
}

/// The oracle's `electron_configuration`, transcribed over the derived tables.
fn configuration(
    order: &[(u32, OrbitalType)],
    exceptions: &[(u32, [Op; 2])],
    z: u32,
) -> Option<Vec<OrbitalFilling>> {
    if z == 0 || z > 118 {
        return None;
    }
    let mut config = madelung(order, z);
    if let Some((_, ops)) = exceptions.iter().find(|x| x.0 == z) {
        for &op in ops {
            set_subshell(&mut config, op);
        }
    }
    config.retain(|e| e.electrons > 0);
    Some(config)
}

/// Each element whose ground state is not the Madelung filling, with the two subshell settings
/// that turn one into the other: changed and emptied subshells in filling order, then a new one.
fn derive_exceptions(order: &[(u32, OrbitalType)]) -> Result<Vec<(u32, [Op; 2])>, String> {
    let mut out = Vec::new();
    for z in 1..=118 {
        let mad = madelung(order, z);
        let want = at::electron_configuration(z).map_err(|e| format!("Z={z}: {e:?}"))?;
        if mad == want {
            continue;
        }
        let same = |a: &OrbitalFilling, b: &OrbitalFilling| a.n == b.n && a.orbital == b.orbital;
        let mut ops: Vec<Op> = Vec::new();
        for x in &mad {
            match want.iter().find(|w| same(w, x)) {
                Some(w) if w.electrons != x.electrons => ops.push((x.n, x.orbital, w.electrons)),
                Some(_) => {}
                None => ops.push((x.n, x.orbital, 0)),
            }
        }
        for w in &want {
            if !mad.iter().any(|x| same(x, w)) {
                ops.push((w.n, w.orbital, w.electrons));
            }
        }
        match ops.as_slice() {
            [a, b] => out.push((z, [*a, *b])),
            _ => return Err(format!("Z={z}: {} subshell changes, not 2", ops.len())),
        }
    }
    Ok(out)
}

/// The oracle's `lamb_shift_nlj_ev`, transcribed over the extracted Bethe logarithms.
#[allow(clippy::cast_precision_loss)]
fn lamb(bethe: &[Arm], z: u32, n: u32, l: u32, two_j: u32) -> Option<f64> {
    if n == 0 || l >= n || (two_j != 2 * l + 1 && two_j + 1 != 2 * l) {
        return None;
    }
    let ln_k0 = bethe.iter().find(|x| x.0 == (n, l))?.1;
    let za = z as f64 * FINE_STRUCTURE;
    let nf = n as f64;
    let pi = core::f64::consts::PI;
    let prefactor = FINE_STRUCTURE / pi * za.powi(4) * ELECTRON_MASS_MEV * 1e6 / (nf * nf * nf);
    let coefficient = if l == 0 {
        let a41 = 4.0 / 3.0 * libm::log(1.0 / (za * za));
        let a40 = -4.0 / 3.0 * ln_k0 + 10.0 / 9.0 - 4.0 / 15.0;
        let a50 = (139.0 / 32.0 - 2.0 * core::f64::consts::LN_2) * pi + 5.0 / 48.0 * pi;
        a41 + a40 + za * a50
    } else {
        let kappa = if two_j == 2 * l + 1 {
            -(l as f64 + 1.0)
        } else {
            l as f64
        };
        -4.0 / 3.0 * ln_k0 - 1.0 / (2.0 * kappa * (2.0 * l as f64 + 1.0))
    };
    Some(prefactor * coefficient)
}

fn member(e: &str, v: &str) -> String {
    format!("TM_{}_{}", snake(e), snake(v))
}

fn letter(o: OrbitalType) -> char {
    o.symbol()
}

pub fn emit(src: &str, provenance: &str) -> Result<String, String> {
    // ── extract / derive, then prove ──
    let order: Vec<(u32, OrbitalType)> =
        u32_variant_pairs(src, "FILLING_ORDER", "OrbitalType", 20)?
            .into_iter()
            .map(|(n, v)| orbital(&v).map(|o| (n, o)))
            .collect::<Result<_, _>>()?;
    let exceptions = derive_exceptions(&order)?;
    if exceptions.len() != 20 {
        return Err(format!(
            "{} configuration exceptions, the oracle documents 20",
            exceptions.len()
        ));
    }
    for z in (0..=120).chain([u32::MAX]) {
        let got = configuration(&order, &exceptions, z);
        if got != at::electron_configuration(z).ok() {
            return Err(format!("configuration proof failed at Z={z}"));
        }
    }
    let bethe = match_f64_arms(src, "bethe_log")?;
    let want_keys: Vec<(u32, u32)> = (1..=4).flat_map(|n| (0..n).map(move |l| (n, l))).collect();
    if bethe.iter().map(|x| x.0).collect::<Vec<_>>() != want_keys {
        return Err("bethe_log: the arms are not (n, l) for n = 1..=4, l < n, in order".into());
    }
    for z in (0..=118).chain([u32::MAX]) {
        for n in 0..=6u32 {
            for l in 0..=n {
                for tj in 0..=2 * n + 1 {
                    let got = lamb(&bethe, z, n, l, tj).map(f64::to_bits);
                    let want = at::lamb_shift_nlj_ev(z, n, l, tj).ok().map(f64::to_bits);
                    if got != want {
                        return Err(format!(
                            "Bethe-logarithm proof failed at Z={z} n={n} l={l} 2j={tj}"
                        ));
                    }
                }
            }
        }
    }
    let ie_cite = z_comments(src, "IONIZATION_ENERGIES_EV", 118)?;
    let ea_cite = z_comments(src, "ELECTRON_AFFINITY_EV", 118)?;
    let symbol = |z: u32| -> String {
        let c = &ie_cite[(z - 1) as usize];
        c.split(" (Z=").next().unwrap_or(c).to_string()
    };

    let mut out = header(
        "atomic data: orbital types, filling order, exceptions, ionization energies, affinities, Bethe logs",
        provenance,
        &[
            "The tables come from the oracle's public API, except three private ones. The filling order and\n\
             the Bethe logarithms are extracted from rust-old/src/atomic.rs; the 20 NIST exceptions are derived\n\
             from electron_configuration. All three are proven against the oracle before this file is written\n\
             (scripts/oracle/src/t_atomic.rs). tests/golden/atomic.tsv pins every use.",
        ],
    );

    // ── enums ──
    let orb_names: Vec<String> = ORBITALS.iter().map(|v| format!("{v:?}")).collect();
    for (rust, variants) in [
        ("OrbitalType", orb_names.clone()),
        (
            "TransitionType",
            vec!["ElectricDipole".to_string(), "Forbidden".to_string()],
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
    out.push_str(&doc_comment(src, "pub enum ElectronAffinity {")?);
    out.push_str(
        "#\n# The kind only: Rust's `Bound(f64)` payload is electron_affinity_ev (port plan §5).\n",
    );
    out.push_str("enum TmElectronAffinity {\n");
    for (k, v) in ["Bound", "Unbound", "Unknown"].iter().enumerate() {
        let _ = writeln!(out, "    {} = {k};", member("ElectronAffinity", v));
    }
    out.push_str("}\nconst TM_ELECTRON_AFFINITY_COUNT = 3;\n\n");

    // ── OrbitalType properties ──
    let ls = ORBITALS
        .iter()
        .map(|o| u8::try_from(o.l()).unwrap_or(u8::MAX))
        .collect();
    out.push_str(&table(
        "TM_ORBITAL_TYPE_L",
        "Orbital angular momentum l. Indexed by TmOrbitalType; mirrors `OrbitalType::l`.",
        &Elem::U8(ls),
        &orb_names,
    )?);
    let me = ORBITALS
        .iter()
        .map(|o| u8::try_from(o.max_electrons()).unwrap_or(u8::MAX))
        .collect();
    out.push_str(&table(
        "TM_ORBITAL_TYPE_MAX_ELECTRONS",
        "Subshell capacity 2(2l + 1). Indexed by TmOrbitalType; mirrors `OrbitalType::max_electrons`.",
        &Elem::U8(me),
        &orb_names,
    )?);
    let sy = ORBITALS
        .iter()
        .map(|o| u8::try_from(u32::from(o.symbol())).unwrap_or(0))
        .collect();
    let sy_lab: Vec<String> = ORBITALS
        .iter()
        .map(|o| format!("'{}'", o.symbol()))
        .collect();
    out.push_str(&table(
        "TM_ORBITAL_TYPE_SYMBOL",
        "The letter symbol, as an ASCII code. Indexed by TmOrbitalType; mirrors `OrbitalType::symbol`.",
        &Elem::U8(sy),
        &sy_lab,
    )?);

    // ── the filling order ──
    out.push_str(&doc_comment(src, "const FILLING_ORDER:")?);
    let _ = writeln!(out, "const _TM_FILLING_ORDER_COUNT = {};\n", order.len());
    let lab: Vec<String> = order
        .iter()
        .map(|&(n, o)| format!("{n}{}", letter(o)))
        .collect();
    let ns = order
        .iter()
        .map(|x| u8::try_from(x.0).unwrap_or(u8::MAX))
        .collect();
    out.push_str(&table(
        "_TM_FILLING_ORDER_N",
        "Column: n.",
        &Elem::U8(ns),
        &lab,
    )?);
    let os = order
        .iter()
        .map(|x| member("OrbitalType", &format!("{:?}", x.1)))
        .collect();
    out.push_str(&table(
        "_TM_FILLING_ORDER_ORBITAL",
        "Column: the orbital type.",
        &Elem::Names(os),
        &lab,
    )?);

    // ── the exceptions ──
    out.push_str(&doc_comment(src, "fn apply_exceptions(")?);
    out.push_str(
        "#\n# Derived from electron_configuration: each element's two subshell settings, applied to the\n\
         # Madelung filling as Rust's set_subshell does (an absent subshell is appended), after\n\
         # which empty subshells are dropped.\n",
    );
    let _ = writeln!(
        out,
        "const _TM_CONFIG_EXCEPTION_COUNT = {};\n",
        exceptions.len()
    );
    let lab: Vec<String> = exceptions
        .iter()
        .map(|(z, ops)| {
            let o: Vec<String> = ops
                .iter()
                .map(|&(n, o, e)| format!("{n}{} → {e}", letter(o)))
                .collect();
            format!("{} (Z={z}): {}", symbol(*z), o.join(", "))
        })
        .collect();
    let zs = exceptions
        .iter()
        .map(|x| u8::try_from(x.0).unwrap_or(u8::MAX))
        .collect();
    out.push_str(&table(
        "_TM_CONFIG_EXCEPTION_Z",
        "Column: Z.",
        &Elem::U8(zs),
        &lab,
    )?);
    for (k, which) in ["1", "2"].iter().enumerate() {
        let n = exceptions
            .iter()
            .map(|x| u8::try_from(x.1[k].0).unwrap_or(u8::MAX))
            .collect();
        out.push_str(&table(
            &format!("_TM_CONFIG_EXCEPTION_N{which}"),
            &format!("Column: setting {which}, n."),
            &Elem::U8(n),
            &lab,
        )?);
        let o = exceptions
            .iter()
            .map(|x| member("OrbitalType", &format!("{:?}", x.1[k].1)))
            .collect();
        out.push_str(&table(
            &format!("_TM_CONFIG_EXCEPTION_ORBITAL{which}"),
            &format!("Column: setting {which}, the orbital type."),
            &Elem::Names(o),
            &lab,
        )?);
        let e = exceptions
            .iter()
            .map(|x| u8::try_from(x.1[k].2).unwrap_or(u8::MAX))
            .collect();
        out.push_str(&table(
            &format!("_TM_CONFIG_EXCEPTION_ELECTRONS{which}"),
            &format!("Column: setting {which}, the electron count."),
            &Elem::U8(e),
            &lab,
        )?);
    }

    // ── ionization energies ──
    let ie = (1..=118)
        .map(|z| at::ionization_energy_ev(z).map_err(|e| format!("Z={z}: {e:?}")))
        .collect::<Result<_, _>>()?;
    out.push_str(&doc_comment(src, "const IONIZATION_ENERGIES_EV:")?);
    out.push_str(&table(
        "_TM_IONIZATION_ENERGY_EV",
        "Indexed by Z − 1. Read through ionization_energy_ev.",
        &Elem::F64(ie),
        &ie_cite,
    )?);

    // ── electron affinities ──
    let mut kinds = Vec::new();
    let mut evs = Vec::new();
    for z in 1..=118 {
        match at::electron_affinity(z).map_err(|e| format!("Z={z}: {e:?}"))? {
            ElectronAffinity::Bound(ev) => {
                kinds.push(member("ElectronAffinity", "Bound"));
                evs.push(ev);
            }
            ElectronAffinity::Unbound => {
                kinds.push(member("ElectronAffinity", "Unbound"));
                evs.push(0.0);
            }
            _ => {
                kinds.push(member("ElectronAffinity", "Unknown"));
                evs.push(0.0);
            }
        }
    }
    out.push_str(&doc_comment(src, "const ELECTRON_AFFINITY_EV:")?);
    out.push_str(&table(
        "_TM_ELECTRON_AFFINITY_KIND",
        "Column: the kind (TmElectronAffinity). Indexed by Z − 1.",
        &Elem::Names(kinds),
        &ea_cite,
    )?);
    out.push_str(&table(
        "_TM_ELECTRON_AFFINITY_EV",
        "Column: the bound value, eV; 0.0 when unbound or unknown. Indexed by Z − 1.",
        &Elem::F64(evs),
        &ea_cite,
    )?);

    // ── Bethe logarithms ──
    out.push_str(&doc_comment(src, "fn bethe_log(")?);
    out.push_str("#\n# Indexed by n(n − 1)/2 + l for n = 1..=4, l < n.\n");
    let lab: Vec<String> = bethe
        .iter()
        .map(|((n, l), _)| format!("{n}{}", ["s", "p", "d", "f"][*l as usize]))
        .collect();
    let vals = bethe.iter().map(|x| x.1).collect();
    out.push_str(&table(
        "_TM_BETHE_LOG",
        "ln k₀(n, l).",
        &Elem::F64(vals),
        &lab,
    )?);
    Ok(out.trim_end().to_string() + "\n")
}
