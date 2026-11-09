//! `src/decay_data.cyr` — the decay modes and the 114 known isotopes (port plan M2).
//!
//! Everything comes from the oracle's API (`known_isotopes()`), and `TmDecayMode` follows the
//! Rust declaration order. Names are NUL-terminated strings in one `u8` byte table with a
//! `u16` offset column. Baked arrays can't hold string pointers, and this way a name is a
//! pointer with no allocation.

use crate::cyr::{Elem, comment_line, doc_comment, header, snake, table};
use crate::m_decay::MODES;
use std::fmt::Write as _;
use tanmatra::decay::known_isotopes;

fn member(v: &str) -> String {
    format!("TM_DECAY_MODE_{}", snake(v))
}

pub fn emit(src: &str, provenance: &str) -> Result<String, String> {
    let mut out = header(
        "decay modes and the known-isotope table (NUBASE2020, 1 y = 365.2422 d)",
        provenance,
        &["Generated from the oracle's known_isotopes(); tests/golden/decay.tsv pins every field."],
    );

    // ── TmDecayMode ──
    out.push_str(&doc_comment(src, "pub enum DecayMode {")?);
    out.push_str("enum TmDecayMode {\n");
    for (k, m) in MODES.iter().enumerate() {
        let _ = writeln!(out, "    {} = {k};", member(&format!("{m:?}")));
    }
    out.push_str("}\n");
    let _ = writeln!(out, "const TM_DECAY_MODE_COUNT = {};\n", MODES.len());

    // ── isotopes ──
    let isos = known_isotopes();
    let lab: Vec<String> = isos.iter().map(|i| i.name.clone()).collect();
    out.push_str(&doc_comment(src, "pub fn known_isotopes()")?);
    let _ = writeln!(out, "const TM_ISOTOPE_COUNT = {};\n", isos.len());
    let z = isos
        .iter()
        .map(|i| u8::try_from(i.nucleus.z()).map_err(|e| e.to_string()))
        .collect::<Result<_, _>>()?;
    out.push_str(&table("TM_ISOTOPE_Z", "Column: Z.", &Elem::U8(z), &lab)?);
    let _ = writeln!(
        out,
        "# Column: A.\nvar TM_ISOTOPE_A: u16[{}] = {{",
        isos.len()
    );
    for i in &isos {
        let _ = writeln!(out, "    {},  # {}", i.nucleus.a(), i.name);
    }
    out.push_str("};\n\n");
    let hl = isos.iter().map(|i| i.half_life_seconds).collect();
    out.push_str(&table(
        "TM_ISOTOPE_HALF_LIFE_SECONDS",
        "Column: half-life, s (+inf: no observed decay).",
        &Elem::F64(hl),
        &lab,
    )?);
    let modes = isos
        .iter()
        .map(|i| member(&format!("{:?}", i.primary_decay)))
        .collect();
    out.push_str(&table(
        "TM_ISOTOPE_PRIMARY_DECAY",
        "Column: the dominant decay branch (TmDecayMode).",
        &Elem::Names(modes),
        &lab,
    )?);
    let iso = isos.iter().map(|i| i.is_isomer).collect();
    out.push_str(&table(
        "TM_ISOTOPE_IS_ISOMER",
        "Column: an isomeric (metastable) state.",
        &Elem::Bool(iso),
        &lab,
    )?);
    let exc = isos.iter().map(|i| i.excitation_energy_kev).collect();
    out.push_str(&table(
        "TM_ISOTOPE_EXCITATION_ENERGY_KEV",
        "Column: excitation energy above the ground state, keV (0 for a ground state).",
        &Elem::F64(exc),
        &lab,
    )?);

    // ── names: offsets + one NUL-terminated byte table ──
    let mut offsets = Vec::new();
    let mut bytes = 0usize;
    for i in &isos {
        offsets.push(bytes);
        bytes += i.name.len() + 1;
    }
    let _ = writeln!(
        out,
        "{}\nvar TM_ISOTOPE_NAME_OFFSET: u16[{}] = {{",
        comment_line("Column: offset of the name in TM_ISOTOPE_NAME_BYTES."),
        isos.len()
    );
    for (o, i) in offsets.iter().zip(&isos) {
        let _ = writeln!(out, "    {o},  # {}", i.name);
    }
    out.push_str("};\n\n");
    let _ = writeln!(
        out,
        "{}\nvar TM_ISOTOPE_NAME_BYTES: u8[{bytes}] = {{",
        comment_line("Every name, NUL-terminated, in table order (ASCII).")
    );
    for i in &isos {
        if !i.name.is_ascii() {
            return Err(format!("isotope name {:?} is not ASCII", i.name));
        }
        let codes: Vec<String> = i.name.bytes().chain([0]).map(|b| b.to_string()).collect();
        let _ = writeln!(out, "    {},  # {}", codes.join(", "), i.name);
    }
    out.push_str("};\n");
    Ok(out)
}
