//! `src/particle_data.cyr` — the Standard Model enums and every tabulated particle
//! property, as baked column arrays indexed by enum value.
//!
//! The enum values follow the Rust declaration order. `e:Quark.Up` in the golden corpus is
//! `TM_QUARK_UP` here: `TM_` + snake(enum) + `_` + snake(variant). Computed properties
//! (`Lepton::decay_width_mev`, `Boson::lifetime_seconds`: ħ over a tabulated value) are not
//! tabulated. `particle.cyr` computes them, and their golden rows are class B.

use crate::cyr::{Elem, doc_comment, doc_comment_after, header, snake, table};
use crate::m_particle::{BOSONS, FORCES, LEPTONS, QUARKS};
use std::fmt::Debug;
use std::fmt::Write as _;

fn names<T: Debug>(vs: &[T]) -> Vec<String> {
    vs.iter().map(|v| format!("{v:?}")).collect()
}

fn member(rust_enum: &str, variant: &str) -> String {
    format!("TM_{}_{}", snake(rust_enum), snake(variant))
}

/// `enum TmX { TM_X_A = 0; .. }` and `const TM_X_COUNT`, under the Rust enum's docs.
fn enum_block(src: &str, rust_enum: &str, variants: &[String]) -> Result<String, String> {
    let mut out = doc_comment(src, &format!("pub enum {rust_enum} {{"))?;
    let _ = writeln!(out, "enum Tm{rust_enum} {{");
    for (k, v) in variants.iter().enumerate() {
        let _ = writeln!(out, "    {} = {k};", member(rust_enum, v));
    }
    out.push_str("}\n");
    let _ = writeln!(
        out,
        "const TM_{}_COUNT = {};\n",
        snake(rust_enum),
        variants.len()
    );
    Ok(out)
}

/// The Rust method's doc line(s), plus which method the table mirrors.
fn method_doc(src: &str, rust_enum: &str, method: &str) -> Result<String, String> {
    let anchor = format!("impl {rust_enum} {{");
    let docs = doc_comment_after(src, &anchor, &format!("pub const fn {method}("))
        .or_else(|_| doc_comment_after(src, &anchor, &format!("pub fn {method}(")))?;
    let first = docs
        .lines()
        .next()
        .unwrap_or("#")
        .trim_start_matches('#')
        .trim();
    Ok(format!(
        "{first}\nIndexed by Tm{rust_enum}; mirrors `{rust_enum}::{method}`."
    ))
}

pub fn emit(src: &str, provenance: &str) -> Result<String, String> {
    let mut out = header(
        "Standard Model particle data (PDG 2024), indexed by enum value",
        provenance,
        &[
            "Enum values follow the Rust declaration order and are append-only. Values are the\n\
           oracle's f64s as round-trip decimals; tests/golden/particle.tsv pins every bit.",
        ],
    );

    let q = names(&QUARKS);
    out.push_str(&enum_block(src, "Quark", &q)?);
    let doc = method_doc(src, "Quark", "mass_mev")?;
    let v = QUARKS.iter().map(|x| x.mass_mev()).collect();
    out.push_str(&table("TM_QUARK_MASS_MEV", &doc, &Elem::F64(v), &q)?);
    let doc = method_doc(src, "Quark", "charge")?;
    let v = QUARKS.iter().map(|x| x.charge()).collect();
    out.push_str(&table("TM_QUARK_CHARGE", &doc, &Elem::F64(v), &q)?);
    let doc = method_doc(src, "Quark", "generation")?;
    let v = QUARKS.iter().map(|x| x.generation()).collect();
    out.push_str(&table("TM_QUARK_GENERATION", &doc, &Elem::U8(v), &q)?);

    let l = names(&LEPTONS);
    out.push_str(&enum_block(src, "Lepton", &l)?);
    let doc = method_doc(src, "Lepton", "mass_mev")?;
    let v = LEPTONS.iter().map(|x| x.mass_mev()).collect();
    out.push_str(&table("TM_LEPTON_MASS_MEV", &doc, &Elem::F64(v), &l)?);
    let doc = method_doc(src, "Lepton", "charge")?;
    let v = LEPTONS.iter().map(|x| x.charge()).collect();
    out.push_str(&table("TM_LEPTON_CHARGE", &doc, &Elem::F64(v), &l)?);
    let doc = method_doc(src, "Lepton", "generation")?;
    let v = LEPTONS.iter().map(|x| x.generation()).collect();
    out.push_str(&table("TM_LEPTON_GENERATION", &doc, &Elem::U8(v), &l)?);
    let doc = method_doc(src, "Lepton", "lifetime_seconds")?;
    let v = LEPTONS.iter().map(|x| x.lifetime_seconds()).collect();
    out.push_str(&table(
        "TM_LEPTON_LIFETIME_SECONDS",
        &doc,
        &Elem::F64(v),
        &l,
    )?);
    let doc = method_doc(src, "Lepton", "is_stable")?;
    let v = LEPTONS.iter().map(|x| x.is_stable()).collect();
    out.push_str(&table("TM_LEPTON_IS_STABLE", &doc, &Elem::Bool(v), &l)?);

    let b = names(&BOSONS);
    out.push_str(&enum_block(src, "Boson", &b)?);
    let doc = method_doc(src, "Boson", "mass_mev")?;
    let v = BOSONS.iter().map(|x| x.mass_mev()).collect();
    out.push_str(&table("TM_BOSON_MASS_MEV", &doc, &Elem::F64(v), &b)?);
    let doc = method_doc(src, "Boson", "spin")?;
    let v = BOSONS.iter().map(|x| x.spin()).collect();
    out.push_str(&table("TM_BOSON_SPIN", &doc, &Elem::U8(v), &b)?);
    let doc = method_doc(src, "Boson", "charge")?;
    let v = BOSONS.iter().map(|x| x.charge()).collect();
    out.push_str(&table("TM_BOSON_CHARGE", &doc, &Elem::F64(v), &b)?);
    let doc = method_doc(src, "Boson", "decay_width_mev")?;
    let v = BOSONS.iter().map(|x| x.decay_width_mev()).collect();
    out.push_str(&table("TM_BOSON_DECAY_WIDTH_MEV", &doc, &Elem::F64(v), &b)?);
    let doc = method_doc(src, "Boson", "is_stable")?;
    let v = BOSONS.iter().map(|x| x.is_stable()).collect();
    out.push_str(&table("TM_BOSON_IS_STABLE", &doc, &Elem::Bool(v), &b)?);

    let fo = names(&FORCES);
    out.push_str(&enum_block(src, "FundamentalForce", &fo)?);
    let doc = method_doc(src, "FundamentalForce", "relative_strength")?;
    let v = FORCES.iter().map(|x| x.relative_strength()).collect();
    out.push_str(&table(
        "TM_FUNDAMENTAL_FORCE_RELATIVE_STRENGTH",
        &doc,
        &Elem::F64(v),
        &fo,
    )?);
    let doc = method_doc(src, "FundamentalForce", "range_meters")?;
    let v = FORCES.iter().map(|x| x.range_meters()).collect();
    out.push_str(&table(
        "TM_FUNDAMENTAL_FORCE_RANGE_METERS",
        &doc,
        &Elem::F64(v),
        &fo,
    )?);

    // mediator() returns a slice per force: flatten it into start / count / members.
    let mut start = Vec::new();
    let mut count = Vec::new();
    let mut members = Vec::new();
    let mut member_labels = Vec::new();
    for force in FORCES {
        let med = force.mediator();
        start.push(u8::try_from(members.len()).map_err(|e| e.to_string())?);
        count.push(u8::try_from(med.len()).map_err(|e| e.to_string())?);
        for boson in med {
            members.push(member("Boson", &format!("{boson:?}")));
            member_labels.push(format!("{force:?}"));
        }
    }
    let doc = method_doc(src, "FundamentalForce", "mediator")?;
    let first = doc.lines().next().unwrap_or_default();
    let flat = format!(
        "{first}\nForce f's mediators are TM_FUNDAMENTAL_FORCE_MEDIATORS[START[f] .. START[f] + COUNT[f]]\n\
         (TmBoson values); mirrors `FundamentalForce::mediator`."
    );
    let doc = format!(
        "{flat}\nStart index into TM_FUNDAMENTAL_FORCE_MEDIATORS, indexed by TmFundamentalForce."
    );
    out.push_str(&table(
        "TM_FUNDAMENTAL_FORCE_MEDIATOR_START",
        &doc,
        &Elem::U8(start),
        &fo,
    )?);
    let doc = "Number of mediators, indexed by TmFundamentalForce.";
    out.push_str(&table(
        "TM_FUNDAMENTAL_FORCE_MEDIATOR_COUNT",
        doc,
        &Elem::U8(count),
        &fo,
    )?);
    let doc = "Every force's mediators, in force order (TmBoson values).";
    out.push_str(&table(
        "TM_FUNDAMENTAL_FORCE_MEDIATORS",
        doc,
        &Elem::Names(members),
        &member_labels,
    )?);

    Ok(out.trim_end().to_string() + "\n")
}
