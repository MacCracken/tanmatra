//! `src/reaction_data.cyr` — the reaction and scattering tables (port plan M4).
//!
//! - **Through the oracle's API:** the seven preset reactions, the thermal cross sections, the
//!   U-235 and Pu-239 fission yields, the two nucleosynthesis pathways, `NucleosynthesisProcess`
//!   and `CNO_NEUTRINO_LOSS_MEV`.
//! - **Extracted, then proven:**
//!   - the private `RESONANCE_INTEGRALS`, by recomputing all three lookups for every (Z, A) with
//!     Z ≤ 130 and A ≤ 300, plus u32::MAX;
//!   - scattering's private `ZETA3`, by recomputing `pair_production_cross_section` on its k ≥ 4
//!     branch for Z = 1..=100 over 13 decades of energy.
//!
//! Strings (preset and pathway names, step reactions) are NUL-terminated UTF-8 in one byte
//! table with a `u16` offset column per string kind.

use crate::cyr::{Elem, comment_line, doc_comment, f64_lit, header, snake, table};
use crate::m_reaction::{PATHWAYS, PRESETS};
use crate::rust_src::{f64_const, f64_tuples};
use std::fmt::Write as _;
use tanmatra::constants::{CLASSICAL_ELECTRON_RADIUS_FM, ELECTRON_MASS_MEV, FINE_STRUCTURE};
use tanmatra::nucleus::Nucleus;
use tanmatra::particle::Lepton;
use tanmatra::reaction::{self as rx, NucleosynthesisProcess};
use tanmatra::scattering as sc;

const PROCESSES: [NucleosynthesisProcess; 3] = [
    NucleosynthesisProcess::SProcess,
    NucleosynthesisProcess::RProcess,
    NucleosynthesisProcess::RpProcess,
];

fn word(n: &Nucleus) -> i64 {
    i64::from(n.z()) | (i64::from(n.a()) << 32)
}

fn lepton_member(l: Lepton) -> String {
    format!("TM_LEPTON_{}", snake(&format!("{l:?}")))
}

/// The scattering oracle's `pair_production_cross_section` on k ≥ 4, transcribed over ζ(3).
fn pair_k4(z: u32, e: f64, zeta3: f64) -> f64 {
    let k = e / ELECTRON_MASS_MEV;
    let pi = core::f64::consts::PI;
    let z2 = f64::from(z) * f64::from(z);
    let prefactor =
        z2 * FINE_STRUCTURE * CLASSICAL_ELECTRON_RADIUS_FM * CLASSICAL_ELECTRON_RADIUS_FM;
    let l = libm::log(2.0 * k);
    let r2 = (2.0 / k) * (2.0 / k);
    let bracket = 28.0 / 9.0 * l - 218.0 / 27.0
        + r2 * (6.0 * l - 3.5 + 2.0 / 3.0 * l * l * l - l * l - pi * pi / 3.0 * l
            + 2.0 * zeta3
            + pi * pi / 6.0)
        - r2 * r2 * (3.0 / 16.0 * l + 0.125)
        - r2 * r2 * r2 * (29.0 / 2304.0 * l - 77.0 / 13_824.0);
    prefactor * bracket
}

/// Appends strings to the shared byte table; returns each one's offset.
fn intern(bytes: &mut Vec<(String, Vec<u8>)>, total: &mut usize, s: &str) -> usize {
    let at = *total;
    let mut b = s.as_bytes().to_vec();
    b.push(0);
    *total += b.len();
    bytes.push((s.to_string(), b));
    at
}

fn offsets_table(name: &str, doc: &str, offs: &[usize], labels: &[String]) -> String {
    let mut out = String::new();
    let _ = writeln!(
        out,
        "{}\nvar {name}: u16[{}] = {{",
        comment_line(doc),
        offs.len()
    );
    for (o, l) in offs.iter().zip(labels) {
        let _ = writeln!(out, "    {o},  # {l}");
    }
    out.push_str("};\n\n");
    out
}

pub fn emit(reaction_src: &str, scattering_src: &str, provenance: &str) -> Result<String, String> {
    // ── extract, then prove ──
    let ri = f64_tuples(reaction_src, "RESONANCE_INTEGRALS", 5)?;
    let mut ri_rows: Vec<(u32, u32, f64, f64, f64)> = Vec::new();
    for t in &ri {
        let z = t[0];
        let a = t[1];
        if z.fract() != 0.0
            || a.fract() != 0.0
            || !(0.0..=255.0).contains(&z)
            || !(0.0..=65535.0).contains(&a)
        {
            return Err(format!(
                "RESONANCE_INTEGRALS: ({z}, {a}) is not a (Z, A) pair"
            ));
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        ri_rows.push((z as u32, a as u32, t[2], t[3], t[4]));
    }
    let look = |z: u32, a: u32, col: usize| -> Option<u64> {
        ri_rows
            .iter()
            .find(|r| r.0 == z && r.1 == a)
            .map(|r| [r.2, r.3, r.4][col].to_bits())
    };
    for z in (0..=130u32).chain([u32::MAX]) {
        for a in (0..=300u32).chain([u32::MAX]) {
            let want = [
                rx::resonance_integral_barns(z, a),
                rx::fission_resonance_integral_barns(z, a),
                rx::n_alpha_resonance_integral_barns(z, a),
            ];
            for (col, w) in want.iter().enumerate() {
                if look(z, a, col) != w.map(f64::to_bits) {
                    return Err(format!(
                        "resonance-integral proof failed at ({z}, {a}), column {col}"
                    ));
                }
            }
        }
    }
    let zeta3 = f64_const(scattering_src, "ZETA3")?;
    for z in 1..=100u32 {
        for k in 0..=130 {
            let e = 4.0 * ELECTRON_MASS_MEV * libm::pow(10.0, f64::from(k) / 10.0);
            if e / ELECTRON_MASS_MEV < 4.0 {
                continue;
            }
            if pair_k4(z, e, zeta3).to_bits() != sc::pair_production_cross_section(z, e).to_bits() {
                return Err(format!("ZETA3 proof failed at Z={z}, E={e}"));
            }
        }
    }

    let mut out = header(
        "reaction and scattering data: presets, thermal cross sections, resonance integrals, yields, pathways",
        provenance,
        &[
            "From the oracle's public API, except two private values, extracted from rust-old/src and proven\n\
             against the oracle before this file is written: the resonance-integral table and scattering's\n\
             ζ(3) (scripts/oracle/src/t_reaction.rs). tests/golden/reaction.tsv and scattering.tsv pin every use.",
        ],
    );

    // ── NucleosynthesisProcess ──
    out.push_str(&doc_comment(
        reaction_src,
        "pub enum NucleosynthesisProcess {",
    )?);
    out.push_str("enum TmNucleosynthesisProcess {\n");
    for (k, p) in PROCESSES.iter().enumerate() {
        let _ = writeln!(
            out,
            "    TM_NUCLEOSYNTHESIS_PROCESS_{} = {k};",
            snake(&format!("{p:?}"))
        );
    }
    out.push_str("}\nconst TM_NUCLEOSYNTHESIS_PROCESS_COUNT = 3;\n\n");

    out.push_str(&doc_comment(
        reaction_src,
        "pub const CNO_NEUTRINO_LOSS_MEV:",
    )?);
    let _ = writeln!(
        out,
        "const TM_CNO_NEUTRINO_LOSS_MEV = {};\n",
        f64_lit(rx::CNO_NEUTRINO_LOSS_MEV)?
    );

    let mut bytes: Vec<(String, Vec<u8>)> = Vec::new();
    let mut total = 0usize;

    // ── presets ──
    let presets: Vec<_> = PRESETS.iter().map(|(n, c)| (*n, c())).collect();
    let lab: Vec<String> = presets.iter().map(|(n, _)| (*n).to_string()).collect();
    out.push_str(
        "# The seven preset reactions (Q from AME2020 atomic masses), indexed by _TM_PRESET_*.\n",
    );
    let _ = writeln!(out, "const _TM_PRESET_COUNT = {};", presets.len());
    for (k, (n, _)) in presets.iter().enumerate() {
        let _ = writeln!(out, "const _TM_PRESET_{} = {k};", n.to_uppercase());
    }
    out.push('\n');
    let offs: Vec<usize> = presets
        .iter()
        .map(|(_, r)| intern(&mut bytes, &mut total, &r.name))
        .collect();
    out.push_str(&offsets_table(
        "_TM_PRESET_NAME_OFFSET",
        "Column: offset of the name in _TM_REACTION_TEXT.",
        &offs,
        &lab,
    ));
    let proj = presets
        .iter()
        .map(|(_, r)| r.projectile.as_ref().map_or(0, word))
        .collect();
    out.push_str(&table(
        "_TM_PRESET_PROJECTILE",
        "Column: the projectile's nucleus word, 0 for none (a neutron-induced reaction).",
        &Elem::I64(proj),
        &lab,
    )?);
    let tgt = presets.iter().map(|(_, r)| word(&r.target)).collect();
    out.push_str(&table(
        "_TM_PRESET_TARGET",
        "Column: the target's nucleus word.",
        &Elem::I64(tgt),
        &lab,
    )?);
    let mut prod_start = Vec::new();
    let mut prod_count = Vec::new();
    let mut products: Vec<i64> = Vec::new();
    let mut prod_lab = Vec::new();
    let mut lep_start = Vec::new();
    let mut lep_count = Vec::new();
    let mut leptons: Vec<String> = Vec::new();
    let mut lep_lab = Vec::new();
    for (n, r) in &presets {
        prod_start.push(u8::try_from(products.len()).map_err(|e| e.to_string())?);
        prod_count.push(u8::try_from(r.products.len()).map_err(|e| e.to_string())?);
        for p in &r.products {
            products.push(word(p));
            prod_lab.push(format!("{n}: Z={} A={}", p.z(), p.a()));
        }
        lep_start.push(u8::try_from(leptons.len()).map_err(|e| e.to_string())?);
        lep_count.push(u8::try_from(r.leptons_out.len()).map_err(|e| e.to_string())?);
        for l in &r.leptons_out {
            leptons.push(lepton_member(*l));
            lep_lab.push((*n).to_string());
        }
    }
    out.push_str(&table(
        "_TM_PRESET_PRODUCT_START",
        "Column: the first product in _TM_REACTION_PRODUCTS.",
        &Elem::U8(prod_start),
        &lab,
    )?);
    out.push_str(&table(
        "_TM_PRESET_PRODUCT_COUNT",
        "Column: the number of product nuclei.",
        &Elem::U8(prod_count),
        &lab,
    )?);
    let nin = presets
        .iter()
        .map(|(_, r)| u8::try_from(r.neutrons_in).unwrap_or(u8::MAX))
        .collect();
    out.push_str(&table(
        "_TM_PRESET_NEUTRONS_IN",
        "Column: free neutrons in the entrance channel.",
        &Elem::U8(nin),
        &lab,
    )?);
    let nout = presets
        .iter()
        .map(|(_, r)| u8::try_from(r.neutrons_out).unwrap_or(u8::MAX))
        .collect();
    out.push_str(&table(
        "_TM_PRESET_NEUTRONS_OUT",
        "Column: free neutrons in the exit channel.",
        &Elem::U8(nout),
        &lab,
    )?);
    out.push_str(&table(
        "_TM_PRESET_LEPTON_START",
        "Column: the first lepton in _TM_REACTION_LEPTONS.",
        &Elem::U8(lep_start),
        &lab,
    )?);
    out.push_str(&table(
        "_TM_PRESET_LEPTON_COUNT",
        "Column: the number of leptons out.",
        &Elem::U8(lep_count),
        &lab,
    )?);
    let q = presets.iter().map(|(_, r)| r.q_value_mev).collect();
    out.push_str(&table(
        "_TM_PRESET_Q_VALUE_MEV",
        "Column: Q, MeV.",
        &Elem::F64(q),
        &lab,
    )?);
    out.push_str(&table(
        "_TM_REACTION_PRODUCTS",
        "Every preset's product nuclei, as nucleus words, in preset order.",
        &Elem::I64(products),
        &prod_lab,
    )?);
    out.push_str(&table(
        "_TM_REACTION_LEPTONS",
        "Every preset's exit-channel leptons (TmLepton; a positron is TM_LEPTON_ELECTRON, as in Rust).",
        &Elem::NamesI64(leptons),
        &lep_lab,
    )?);

    // ── thermal cross sections ──
    let th = rx::thermal_neutron_cross_sections();
    let lab: Vec<String> = th
        .iter()
        .map(|x| format!("Z={} A={}", x.target.z(), x.target.a()))
        .collect();
    out.push_str(&doc_comment(
        reaction_src,
        "pub fn thermal_neutron_cross_sections()",
    )?);
    let _ = writeln!(out, "const TM_THERMAL_XS_COUNT = {};\n", th.len());
    let zs = th
        .iter()
        .map(|x| u8::try_from(x.target.z()).unwrap_or(u8::MAX))
        .collect();
    out.push_str(&table(
        "TM_THERMAL_XS_Z",
        "Column: target Z.",
        &Elem::U8(zs),
        &lab,
    )?);
    let as_ = th
        .iter()
        .map(|x| u8::try_from(x.target.a()).unwrap_or(u8::MAX))
        .collect();
    out.push_str(&table(
        "TM_THERMAL_XS_A",
        "Column: target A.",
        &Elem::U8(as_),
        &lab,
    )?);
    for (col, doc, get) in [
        (
            "ABSORPTION_BARNS",
            "Column: absorption, b (capture + fission + (n,α) + (n,p)).",
            0,
        ),
        ("FISSION_BARNS", "Column: fission, b.", 1),
        (
            "SCATTERING_BARNS",
            "Column: free-atom elastic scattering, b.",
            2,
        ),
    ] {
        let v = th
            .iter()
            .map(|x| [x.absorption_barns, x.fission_barns, x.scattering_barns][get])
            .collect();
        out.push_str(&table(
            &format!("TM_THERMAL_XS_{col}"),
            doc,
            &Elem::F64(v),
            &lab,
        )?);
    }

    // ── resonance integrals ──
    let lab: Vec<String> = ri_rows
        .iter()
        .map(|r| format!("Z={} A={}", r.0, r.1))
        .collect();
    out.push_str(&doc_comment(reaction_src, "const RESONANCE_INTEGRALS:")?);
    let _ = writeln!(
        out,
        "const _TM_RESONANCE_INTEGRAL_COUNT = {};\n",
        ri_rows.len()
    );
    let zs = ri_rows
        .iter()
        .map(|r| u8::try_from(r.0).unwrap_or(u8::MAX))
        .collect();
    out.push_str(&table(
        "_TM_RESONANCE_INTEGRAL_Z",
        "Column: Z.",
        &Elem::U8(zs),
        &lab,
    )?);
    let as_ = ri_rows
        .iter()
        .map(|r| u8::try_from(r.1).unwrap_or(u8::MAX))
        .collect();
    out.push_str(&table(
        "_TM_RESONANCE_INTEGRAL_A",
        "Column: A.",
        &Elem::U8(as_),
        &lab,
    )?);
    for (col, doc, k) in [
        ("CAPTURE", "Column: capture I_γ, b.", 0),
        ("FISSION", "Column: fission I_f, b.", 1),
        ("N_ALPHA", "Column: (n,α) I_α, b.", 2),
    ] {
        let v = ri_rows.iter().map(|r| [r.2, r.3, r.4][k]).collect();
        out.push_str(&table(
            &format!("_TM_RESONANCE_INTEGRAL_{col}"),
            doc,
            &Elem::F64(v),
            &lab,
        )?);
    }

    // ── fission yields ──
    for (tag, fname, ys) in [
        (
            "U235",
            "pub fn u235_fission_yields()",
            rx::u235_fission_yields(),
        ),
        (
            "PU239",
            "pub fn pu239_fission_yields()",
            rx::pu239_fission_yields(),
        ),
    ] {
        let lab: Vec<String> = ys.iter().map(|y| format!("A={}", y.mass_number)).collect();
        out.push_str(&doc_comment(reaction_src, fname)?);
        let _ = writeln!(out, "const TM_{tag}_FISSION_YIELD_COUNT = {};\n", ys.len());
        let ms = ys
            .iter()
            .map(|y| u8::try_from(y.mass_number).unwrap_or(u8::MAX))
            .collect();
        out.push_str(&table(
            &format!("TM_{tag}_FISSION_YIELD_MASS_NUMBER"),
            "Column: mass number A.",
            &Elem::U8(ms),
            &lab,
        )?);
        let fr = ys.iter().map(|y| y.yield_fraction).collect();
        out.push_str(&table(
            &format!("TM_{tag}_FISSION_YIELD_FRACTION"),
            "Column: chain yield, per fission.",
            &Elem::F64(fr),
            &lab,
        )?);
    }

    // ── nucleosynthesis pathways ──
    let paths: Vec<_> = PATHWAYS.iter().map(|(n, c)| (*n, c())).collect();
    let lab: Vec<String> = paths.iter().map(|(n, _)| (*n).to_string()).collect();
    out.push_str("# The two nucleosynthesis pathways, indexed by _TM_PATHWAY_*.\n");
    for (k, (n, _)) in paths.iter().enumerate() {
        let _ = writeln!(out, "const _TM_PATHWAY_{} = {k};", n.to_uppercase());
    }
    out.push('\n');
    let offs: Vec<usize> = paths
        .iter()
        .map(|(_, p)| intern(&mut bytes, &mut total, &p.name))
        .collect();
    out.push_str(&offsets_table(
        "_TM_PATHWAY_NAME_OFFSET",
        "Column: offset of the name in _TM_REACTION_TEXT.",
        &offs,
        &lab,
    ));
    let procs = paths
        .iter()
        .map(|(_, p)| {
            format!(
                "TM_NUCLEOSYNTHESIS_PROCESS_{}",
                snake(&format!("{:?}", p.process))
            )
        })
        .collect();
    out.push_str(&table(
        "_TM_PATHWAY_PROCESS",
        "Column: the process (TmNucleosynthesisProcess).",
        &Elem::Names(procs),
        &lab,
    )?);
    let seeds = paths.iter().map(|(_, p)| word(&p.seed)).collect();
    out.push_str(&table(
        "_TM_PATHWAY_SEED",
        "Column: the seed's nucleus word.",
        &Elem::I64(seeds),
        &lab,
    )?);
    let mut st_start = Vec::new();
    let mut st_count = Vec::new();
    let mut st_nuc = Vec::new();
    let mut st_off = Vec::new();
    let mut st_lab = Vec::new();
    for (n, p) in &paths {
        st_start.push(u8::try_from(st_nuc.len()).map_err(|e| e.to_string())?);
        st_count.push(u8::try_from(p.steps.len()).map_err(|e| e.to_string())?);
        for st in &p.steps {
            st_nuc.push(word(&st.nucleus));
            st_off.push(intern(&mut bytes, &mut total, &st.reaction));
            st_lab.push(format!("{n}: {}", st.reaction));
        }
    }
    out.push_str(&table(
        "_TM_PATHWAY_STEP_START",
        "Column: the first step in _TM_PATHWAY_STEP_*.",
        &Elem::U8(st_start),
        &lab,
    )?);
    out.push_str(&table(
        "_TM_PATHWAY_STEP_COUNT",
        "Column: the number of steps.",
        &Elem::U8(st_count),
        &lab,
    )?);
    out.push_str(&table(
        "_TM_PATHWAY_STEP_NUCLEUS",
        "Every step's nucleus word, in pathway order.",
        &Elem::I64(st_nuc),
        &st_lab,
    )?);
    out.push_str(&offsets_table(
        "_TM_PATHWAY_STEP_REACTION_OFFSET",
        "Every step's reaction text: its offset in _TM_REACTION_TEXT.",
        &st_off,
        &st_lab,
    ));

    // ── the shared text ──
    let _ = writeln!(
        out,
        "{}\nvar _TM_REACTION_TEXT: u8[{total}] = {{",
        comment_line(
            "Every name and step reaction above, NUL-terminated UTF-8, in the order of their offsets."
        )
    );
    for (label, b) in &bytes {
        for (k, chunk) in b.chunks(12).enumerate() {
            let codes: Vec<String> = chunk.iter().map(u8::to_string).collect();
            if k == 0 {
                let _ = writeln!(out, "    {},  # {label}", codes.join(", "));
            } else {
                let _ = writeln!(out, "    {},", codes.join(", "));
            }
        }
    }
    out.push_str("};\n\n");

    // ── scattering's ζ(3) ──
    out.push_str(&doc_comment(scattering_src, "const ZETA3:")?);
    let _ = writeln!(out, "const _TM_ZETA3 = {};", f64_lit(zeta3)?);
    Ok(out.trim_end().to_string() + "\n")
}
