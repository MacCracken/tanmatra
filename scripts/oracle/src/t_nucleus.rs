//! `src/nucleus_data.cyr` — the nuclear data tables and model coefficients (port plan M2).
//!
//! - **Public tables, enumerated through the oracle's API:**
//!   - AME2020 mass excesses, charge radii and moments: the (Z, A) lookups, over Z = 1..120
//!     and A = Z..300;
//!   - the superallowed transitions (`superallowed_ft_values`);
//!   - the Mayer–Jensen levels (`shell_model_levels`);
//!   - the shell closures (`next_shell_closure`);
//!   - the average Ft (`superallowed_average_ft`).
//! - **Private values, through `rust_src` and each proven below** by recomputing what the
//!   oracle computes from it, over the whole M2 grid, to identical bits:
//!   - the SEMF coefficients (proof: `binding_energy`);
//!   - the Myers–Swiatecki C, c and magic numbers (proof: `binding_energy_shell_corrected`);
//!   - AME2020's u→keV factor (proof: `experimental_atomic_mass_amu`);
//!   - the empirical proton and neutron level orders (proof: `ground_state_spin_parity`).

use crate::cyr::{Elem, comment_line, doc_comment, f64_lit, header, table};
use crate::m_nucleus::grid;
use crate::rust_src::{f64_const, u32_array, u32_triples};
use std::fmt::Write as _;
use tanmatra::nucleus::{
    Nucleus, ground_state_spin_parity, next_shell_closure, shell_model_levels,
    superallowed_average_ft, superallowed_ft_values,
};

struct Coeffs {
    a_v: f64,
    a_s: f64,
    a_c: f64,
    a_a: f64,
    a_p: f64,
    ms_c: f64,
    ms_small_c: f64,
    ms_magic: Vec<u32>,
}

/// The oracle's `binding_energy`, transcribed operation for operation (the proof of the
/// SEMF coefficients).
fn semf(k: &Coeffs, z: u32, a: u32) -> f64 {
    if a == 1 {
        return 0.0;
    }
    let (af, zf) = (f64::from(a), f64::from(z));
    let c = libm::cbrt(af);
    let volume = k.a_v * af;
    let surface = k.a_s * (c * c);
    let coulomb = k.a_c * zf * (zf - 1.0) / c;
    let asym = k.a_a * ((af - 2.0 * zf) * (af - 2.0 * zf)) / af;
    let n = a - z;
    let denom = libm::sqrt(af);
    let delta = if denom == 0.0 {
        0.0
    } else if z.is_multiple_of(2) && n.is_multiple_of(2) {
        k.a_p / denom
    } else if !z.is_multiple_of(2) && !n.is_multiple_of(2) {
        -k.a_p / denom
    } else {
        0.0
    };
    volume - surface - coulomb - asym + delta
}

fn five_thirds(x: f64) -> f64 {
    let c = libm::cbrt(x);
    x * c * c
}

fn ms_f(k: &Coeffs, count: u32) -> f64 {
    for w in k.ms_magic.windows(2) {
        let (lo, hi) = (w[0], w[1]);
        if count <= hi {
            let (lo_f, hi_f, n) = (f64::from(lo), f64::from(hi), f64::from(count));
            let lo53 = five_thirds(lo_f);
            let q = 0.6 * (five_thirds(hi_f) - lo53) / (hi_f - lo_f);
            return q * (n - lo_f) - 0.6 * (five_thirds(n) - lo53);
        }
    }
    0.0
}

/// The oracle's `binding_energy_shell_corrected`, transcribed (the proof of C, c, MS_MAGIC).
fn shell_corrected(k: &Coeffs, z: u32, a: u32) -> f64 {
    if a == 1 {
        return 0.0;
    }
    let n = a - z;
    let af = f64::from(a);
    let c = libm::cbrt(af / 2.0);
    let s = k.ms_c * ((ms_f(k, n) + ms_f(k, z)) / (c * c) - k.ms_small_c * libm::cbrt(af));
    semf(k, z, a) - s
}

type Order = [(u32, u32, u32)];

fn last_level(order: &Order, count: u32) -> ((u32, u32, u32), u32) {
    let mut rem = count;
    for &lv in order {
        let deg = lv.2 + 1;
        if rem <= deg {
            return (lv, rem);
        }
        rem -= deg;
    }
    (order[order.len() - 1], 1)
}

/// The oracle's `ground_state_spin_parity`, transcribed (the proof of the level orders).
fn spin_parity(p_ord: &Order, n_ord: &Order, z: u32, n: u32) -> (u32, i32) {
    let par = |l: u32| if l.is_multiple_of(2) { 1 } else { -1 };
    match (!z.is_multiple_of(2), !n.is_multiple_of(2)) {
        (false, false) => (0, 1),
        (true, false) => {
            let (p, _) = last_level(p_ord, z);
            (p.2, par(p.1))
        }
        (false, true) => {
            let (q, _) = last_level(n_ord, n);
            (q.2, par(q.1))
        }
        (true, true) => {
            let (p, po) = last_level(p_ord, z);
            let (q, qo) = last_level(n_ord, n);
            // degeneracy / 2, as the oracle writes it (`p_occ > p.degeneracy() / 2`).
            let (p_deg, q_deg) = (p.2 + 1, q.2 + 1);
            let ph = po > p_deg / 2;
            let qh = qo > q_deg / 2;
            let sp = i64::from(p.2) - 2 * i64::from(p.1);
            let sq = i64::from(q.2) - 2 * i64::from(q.1);
            let tj = if ph != qh {
                p.2 + q.2 - 2
            } else if sp + sq == 0 {
                u32::try_from((i64::from(p.2) - i64::from(q.2)).abs()).unwrap_or(0)
            } else {
                p.2 + q.2
            };
            (tj, par(p.1 + q.1))
        }
    }
}

fn check(ok: bool, what: &str) -> Result<(), String> {
    if ok {
        Ok(())
    } else {
        Err(format!("extraction proof failed: {what}"))
    }
}

/// Every (Z, A) the oracle tabulates for `f`, over Z = 1..120, A = Z..300.
fn enumerate<T>(f: impl Fn(&Nucleus) -> Option<T>) -> Vec<(u32, u32, T)> {
    let mut v = Vec::new();
    for z in 1..=120 {
        for a in z..=300 {
            if let Ok(n) = Nucleus::new(z, a)
                && let Some(x) = f(&n)
            {
                v.push((z, a, x));
            }
        }
    }
    v
}

fn za_labels<T>(rows: &[(u32, u32, T)]) -> Vec<String> {
    rows.iter()
        .map(|(z, a, _)| format!("Z={z} A={a}"))
        .collect()
}

fn u8s(v: impl Iterator<Item = u32>) -> Result<Vec<u8>, String> {
    v.map(|x| u8::try_from(x).map_err(|e| e.to_string()))
        .collect()
}

/// A u16 column: written through `Elem::U8` would truncate, so emit it as text directly.
fn u16_table(name: &str, doc: &str, vals: &[u32], labels: &[String]) -> String {
    let mut out = String::new();
    for line in doc.lines() {
        let _ = writeln!(out, "{}", comment_line(line));
    }
    let _ = writeln!(out, "var {name}: u16[{}] = {{", vals.len());
    for (v, l) in vals.iter().zip(labels) {
        let _ = writeln!(out, "    {v},  # {l}");
    }
    out.push_str("};\n\n");
    out
}

pub fn emit(src: &str, provenance: &str) -> Result<String, String> {
    // ── private values, extracted then proven ──
    let k = Coeffs {
        a_v: f64_const(src, "A_V")?,
        a_s: f64_const(src, "A_S")?,
        a_c: f64_const(src, "A_C")?,
        a_a: f64_const(src, "A_A")?,
        a_p: f64_const(src, "A_P")?,
        ms_c: f64_const(src, "MS_SHELL_C")?,
        ms_small_c: f64_const(src, "MS_SHELL_SMALL_C")?,
        ms_magic: u32_array(src, "MS_MAGIC", 9)?,
    };
    let u_kev = f64_const(src, "AME2020_U_KEV")?;
    let p_ord = u32_triples(src, "PROTON_LEVEL_ORDER", 32)?;
    let n_ord = u32_triples(src, "NEUTRON_LEVEL_ORDER", 32)?;
    for (z, a) in grid() {
        let nuc = Nucleus::new(z, a).map_err(|e| e.to_string())?;
        let be = nuc.binding_energy().to_bits();
        check(
            semf(&k, z, a).to_bits() == be,
            &format!("binding_energy({z}, {a})"),
        )?;
        let sc = nuc.binding_energy_shell_corrected().to_bits();
        check(
            shell_corrected(&k, z, a).to_bits() == sc,
            &format!("shell_corrected({z}, {a})"),
        )?;
        let want = ground_state_spin_parity(&nuc);
        check(
            spin_parity(&p_ord, &n_ord, z, a - z) == want,
            &format!("J^pi({z}, {a})"),
        )?;
    }
    let ame = enumerate(Nucleus::experimental_mass_excess_kev);
    for &(z, a, me) in &ame {
        let nuc = Nucleus::new(z, a).map_err(|e| e.to_string())?;
        let want = nuc.experimental_atomic_mass_amu().map(f64::to_bits);
        check(
            Some((f64::from(a) + me / u_kev).to_bits()) == want,
            &format!("AME u({z}, {a})"),
        )?;
    }

    let mut out = header(
        "nuclear data: AME2020, charge radii, moments, superallowed ft, shell-model levels, model coefficients",
        provenance,
        &[
            "Public tables are enumerated through the oracle's API. The private coefficients and level\n\
           orders are extracted from rust-old/src/nucleus.rs and proven against the oracle over the\n\
           whole M2 grid before this file is written (scripts/oracle/src/t_nucleus.rs).\n\
           tests/golden/nucleus.tsv pins every value's use.",
        ],
    );

    // ── model coefficients (private in Rust → `_TM_`) ──
    let nsrc = src;
    out.push_str(
        "# Bethe–Weizsäcker coefficients (MeV): least-squares fit to the 2484 AME2020 binding\n",
    );
    out.push_str("# energies with A >= 16 (rust-old/src/nucleus.rs A_V .. A_P).\n");
    for (name, v) in [
        ("A_V", k.a_v),
        ("A_S", k.a_s),
        ("A_C", k.a_c),
        ("A_A", k.a_a),
        ("A_P", k.a_p),
    ] {
        let _ = writeln!(out, "const _TM_SEMF_{name} = {};", f64_lit(v)?);
    }
    out.push('\n');
    out.push_str(&doc_comment(nsrc, "const MS_SHELL_C:")?);
    let _ = writeln!(out, "const _TM_MS_SHELL_C = {};\n", f64_lit(k.ms_c)?);
    out.push_str(&doc_comment(nsrc, "const MS_SHELL_SMALL_C:")?);
    let _ = writeln!(
        out,
        "const _TM_MS_SHELL_SMALL_C = {};\n",
        f64_lit(k.ms_small_c)?
    );
    out.push_str(&doc_comment(nsrc, "const AME2020_U_KEV:")?);
    let _ = writeln!(out, "const _TM_AME2020_U_KEV = {};\n", f64_lit(u_kev)?);
    let (avg, unc) = superallowed_average_ft();
    out.push_str("# World-average corrected Ft and its total uncertainty, s (Hardy & Towner 2020, eq. 22).\n");
    let _ = writeln!(
        out,
        "const TM_SUPERALLOWED_AVERAGE_FT_SECONDS = {};",
        f64_lit(avg)?
    );
    let _ = writeln!(
        out,
        "const TM_SUPERALLOWED_AVERAGE_FT_UNCERTAINTY = {};\n",
        f64_lit(unc)?
    );

    // ── AME2020 ──
    let doc = doc_comment(nsrc, "const AME2020_MASS_EXCESS:")?;
    let doc = doc
        .lines()
        .map(|l| l.trim_start_matches('#').trim_start())
        .collect::<Vec<_>>()
        .join("\n");
    let lab = za_labels(&ame);
    let _ = writeln!(out, "const TM_AME2020_COUNT = {};\n", ame.len());
    out.push_str(&table(
        "TM_AME2020_Z",
        &format!("{doc}\nColumn: Z."),
        &Elem::U8(u8s(ame.iter().map(|r| r.0))?),
        &lab,
    )?);
    out.push_str(&u16_table(
        "TM_AME2020_A",
        "Column: A.",
        &ame.iter().map(|r| r.1).collect::<Vec<_>>(),
        &lab,
    ));
    out.push_str(&table(
        "TM_AME2020_MASS_EXCESS_KEV",
        "Column: mass excess, keV.",
        &Elem::F64(ame.iter().map(|r| r.2).collect()),
        &lab,
    )?);

    // ── charge radii ──
    let radii = enumerate(Nucleus::charge_radius_fm);
    let lab = za_labels(&radii);
    let doc = doc_comment(nsrc, "const CHARGE_RADII:")?;
    let doc = doc
        .lines()
        .map(|l| l.trim_start_matches('#').trim_start())
        .collect::<Vec<_>>()
        .join("\n");
    let _ = writeln!(out, "const TM_CHARGE_RADII_COUNT = {};\n", radii.len());
    out.push_str(&table(
        "TM_CHARGE_RADII_Z",
        &format!("{doc}\nColumn: Z."),
        &Elem::U8(u8s(radii.iter().map(|r| r.0))?),
        &lab,
    )?);
    out.push_str(&u16_table(
        "TM_CHARGE_RADII_A",
        "Column: A.",
        &radii.iter().map(|r| r.1).collect::<Vec<_>>(),
        &lab,
    ));
    out.push_str(&table(
        "TM_CHARGE_RADII_FM",
        "Column: RMS charge radius, fm.",
        &Elem::F64(radii.iter().map(|r| r.2).collect()),
        &lab,
    )?);

    // ── moments ──
    let mom = enumerate(Nucleus::nuclear_moments);
    let lab = za_labels(&mom);
    let doc = doc_comment(nsrc, "const NUCLEAR_MOMENTS:")?;
    let doc = doc
        .lines()
        .map(|l| l.trim_start_matches('#').trim_start())
        .collect::<Vec<_>>()
        .join("\n");
    let _ = writeln!(out, "const TM_NUCLEAR_MOMENTS_COUNT = {};\n", mom.len());
    out.push_str(&table(
        "TM_NUCLEAR_MOMENTS_Z",
        &format!("{doc}\nColumn: Z."),
        &Elem::U8(u8s(mom.iter().map(|r| r.0))?),
        &lab,
    )?);
    out.push_str(&u16_table(
        "TM_NUCLEAR_MOMENTS_A",
        "Column: A.",
        &mom.iter().map(|r| r.1).collect::<Vec<_>>(),
        &lab,
    ));
    let mu = mom.iter().map(|r| r.2.magnetic_dipole_mu_n).collect();
    out.push_str(&table(
        "TM_NUCLEAR_MOMENTS_MU_N",
        "Column: magnetic dipole moment, nuclear magnetons.",
        &Elem::F64(mu),
        &lab,
    )?);
    let q = mom.iter().map(|r| r.2.electric_quadrupole_barn).collect();
    out.push_str(&table(
        "TM_NUCLEAR_MOMENTS_Q_BARN",
        "Column: electric quadrupole moment, barns.",
        &Elem::F64(q),
        &lab,
    )?);

    // ── superallowed ──
    let sa = superallowed_ft_values();
    let lab: Vec<String> = sa
        .iter()
        .map(|d| {
            format!(
                "Z={} A={} -> Z={} A={}",
                d.parent.z(),
                d.parent.a(),
                d.daughter.z(),
                d.daughter.a()
            )
        })
        .collect();
    let doc = doc_comment(nsrc, "const SUPERALLOWED_FT_VALUES:")?;
    let doc = doc
        .lines()
        .map(|l| l.trim_start_matches('#').trim_start())
        .collect::<Vec<_>>()
        .join("\n");
    let _ = writeln!(out, "const TM_SUPERALLOWED_COUNT = {};\n", sa.len());
    out.push_str(&table(
        "TM_SUPERALLOWED_PARENT_Z",
        &format!("{doc}\nColumn: parent Z."),
        &Elem::U8(u8s(sa.iter().map(|d| d.parent.z()))?),
        &lab,
    )?);
    out.push_str(&u16_table(
        "TM_SUPERALLOWED_PARENT_A",
        "Column: parent A.",
        &sa.iter().map(|d| d.parent.a()).collect::<Vec<_>>(),
        &lab,
    ));
    out.push_str(&table(
        "TM_SUPERALLOWED_DAUGHTER_Z",
        "Column: daughter Z.",
        &Elem::U8(u8s(sa.iter().map(|d| d.daughter.z()))?),
        &lab,
    )?);
    out.push_str(&u16_table(
        "TM_SUPERALLOWED_DAUGHTER_A",
        "Column: daughter A.",
        &sa.iter().map(|d| d.daughter.a()).collect::<Vec<_>>(),
        &lab,
    ));
    out.push_str(&table(
        "TM_SUPERALLOWED_FT_SECONDS",
        "Column: measured ft, s.",
        &Elem::F64(sa.iter().map(|d| d.ft_seconds).collect()),
        &lab,
    )?);
    let big: Vec<f64> = sa
        .iter()
        .map(|d| {
            d.corrected_ft_seconds()
                .ok_or("a superallowed row has no corrected Ft")
        })
        .collect::<Result<_, _>>()?;
    out.push_str(&table(
        "TM_SUPERALLOWED_CORRECTED_FT_SECONDS",
        "Column: corrected Ft = ft(1 + δ'_R)(1 + δ_NS − δ_C), s.",
        &Elem::F64(big),
        &lab,
    )?);

    // ── shell model ──
    let ms = &k.ms_magic;
    let lab: Vec<String> = (0..ms.len()).map(|i| format!("M_{i}")).collect();
    out.push_str(&u16_table(
        "_TM_MS_MAGIC",
        "Myers–Swiatecki magic numbers M_0 .. M_8 (rust-old MS_MAGIC).",
        ms,
        &lab,
    ));
    let mut closures = Vec::new();
    let mut c = 0;
    while let Some(nx) = next_shell_closure(c) {
        closures.push(nx);
        c = nx;
    }
    let lab: Vec<String> = closures.iter().map(|c| format!("{c}")).collect();
    let _ = writeln!(out, "const TM_SHELL_CLOSURE_COUNT = {};\n", closures.len());
    out.push_str(&u16_table(
        "TM_SHELL_CLOSURES",
        "Shell closures (magic numbers), ascending (next_shell_closure).",
        &closures,
        &lab,
    ));

    let levels = shell_model_levels();
    let lab: Vec<String> = levels
        .iter()
        .map(|&(n, l, j)| format!("{n} l={l} 2j={j}"))
        .collect();
    let doc = doc_comment(nsrc, "const SHELL_MODEL_LEVELS:")?;
    let doc = doc
        .lines()
        .map(|l| l.trim_start_matches('#').trim_start())
        .collect::<Vec<_>>()
        .join("\n");
    let _ = writeln!(out, "const TM_SHELL_LEVEL_COUNT = {};\n", levels.len());
    for (col, idx, d) in [
        ("N_SHELL", 0, format!("{doc}\nColumn: n_shell.")),
        ("L", 1, "Column: l.".into()),
        ("TWO_J", 2, "Column: 2j.".into()),
    ] {
        let v = levels.iter().map(|t| [t.0, t.1, t.2][idx]);
        out.push_str(&table(
            &format!("TM_SHELL_LEVEL_{col}"),
            &d,
            &Elem::U8(u8s(v)?),
            &lab,
        )?);
    }
    for (name, ord, rust) in [
        ("PROTON", &p_ord, "PROTON_LEVEL_ORDER"),
        ("NEUTRON", &n_ord, "NEUTRON_LEVEL_ORDER"),
    ] {
        let lab: Vec<String> = ord
            .iter()
            .map(|&(n, l, j)| format!("{n} l={l} 2j={j}"))
            .collect();
        let doc = doc_comment(nsrc, &format!("const {rust}:"))?;
        let doc = doc
            .lines()
            .map(|l| l.trim_start_matches('#').trim_start())
            .collect::<Vec<_>>()
            .join("\n");
        for (col, idx, d) in [
            ("N_SHELL", 0, format!("{doc}\nColumn: n_shell.")),
            ("L", 1, "Column: l.".into()),
            ("TWO_J", 2, "Column: 2j.".into()),
        ] {
            let v = ord.iter().map(|t| [t.0, t.1, t.2][idx]);
            out.push_str(&table(
                &format!("_TM_{name}_ORDER_{col}"),
                &d,
                &Elem::U8(u8s(v)?),
                &lab,
            )?);
        }
    }
    Ok(out.trim_end().to_string() + "\n")
}
