//! Cyrius source emission for the generated data modules (port plan §2 step 6).
//!
//! Rules the output must keep. Each was probed on the 6.7.18 toolchain.
//! - **f64 values are decimal literals, written out positionally.** The lexer has no
//!   exponent notation. Decimals have been correctly rounded since 6.6.10, so the shortest
//!   round-trip digits give back the oracle's exact bits.
//! - **Hex is never used in a `const`**, where it is an integer. A full-width hex pattern in
//!   a typed `f64` element keeps its bits without a warning, but decimals keep every value
//!   readable in review, so tables use decimals too.
//! - **±∞ are `(1.0 / 0.0)` and `-(1.0 / 0.0)`.** NaN can't be a compile-time f64 (the
//!   compiler refuses it: its bits depend on the target), so a NaN in a table is an error
//!   here.
//! - **Every global gets the `TM_` prefix.** A bundle is concatenated into its consumers, so
//!   names share one namespace.

use std::fmt::Write as _;

/// The Cyrius literal for `x`, checked to parse back to the same bits.
pub fn f64_lit(x: f64) -> Result<String, String> {
    if x.is_nan() {
        return Err("NaN cannot be a compile-time f64 in Cyrius".to_string());
    }
    if x.is_infinite() {
        return Ok(if x > 0.0 {
            "(1.0 / 0.0)"
        } else {
            "-(1.0 / 0.0)"
        }
        .to_string());
    }
    let lit = positional(x);
    let back: f64 = lit.parse().map_err(|e| format!("{lit}: {e}"))?;
    if back.to_bits() != x.to_bits() {
        return Err(format!("{lit} does not round-trip to {:016x}", x.to_bits()));
    }
    Ok(lit)
}

/// The shortest round-trip digits of `x` (Rust's `{:e}`), written without an exponent.
fn positional(x: f64) -> String {
    let sci = format!("{:e}", x.abs());
    let (mant, exp) = sci.split_once('e').unwrap_or((&sci, "0"));
    let exp: i64 = exp.parse().unwrap_or(0);
    let digits: String = mant.chars().filter(char::is_ascii_digit).collect();
    let point = exp + 1; // digits before the decimal point
    let n = i64::try_from(digits.len()).unwrap_or(i64::MAX);
    let body = if point <= 0 {
        let zeros = "0".repeat(usize::try_from(-point).unwrap_or(0));
        format!("0.{zeros}{digits}")
    } else if point >= n {
        let zeros = "0".repeat(usize::try_from(point - n).unwrap_or(0));
        format!("{digits}{zeros}.0")
    } else {
        let at = usize::try_from(point).unwrap_or(0);
        format!("{}.{}", &digits[..at], &digits[at..])
    };
    if x.is_sign_negative() {
        format!("-{body}")
    } else {
        body
    }
}

/// `ElectronNeutrino` → `ELECTRON_NEUTRINO`, `WPlus` → `W_PLUS`, `TAI` → `TAI`,
/// `Cesium133` → `CESIUM133`. An underscore goes before a capital that follows a lowercase
/// letter or digit, or that starts a new capitalised word after a run of capitals. A run of
/// capitals (an acronym) stays one word.
pub fn snake(camel: &str) -> String {
    let c: Vec<char> = camel.chars().collect();
    let mut out = String::new();
    for (k, &ch) in c.iter().enumerate() {
        if ch.is_ascii_uppercase() && k > 0 {
            let prev = c[k - 1];
            let next_lower = c.get(k + 1).is_some_and(char::is_ascii_lowercase);
            if prev.is_ascii_lowercase()
                || prev.is_ascii_digit()
                || (prev.is_ascii_uppercase() && next_lower)
            {
                out.push('_');
            }
        }
        out.push(ch.to_ascii_uppercase());
    }
    out
}

/// The `///` doc block directly above the first line of `src` that starts with `item`
/// (attributes in between are skipped), as `# ` comment lines. Empty if there is none.
/// Only comments are read from the Rust source; values always come from the oracle's API.
pub fn doc_comment(src: &str, item: &str) -> Result<String, String> {
    doc_comment_after(src, "", item)
}

/// As [`doc_comment`], for the first `item` after the first line starting with `anchor`
/// (e.g. a method inside `impl Quark {`).
pub fn doc_comment_after(src: &str, anchor: &str, item: &str) -> Result<String, String> {
    let lines: Vec<&str> = src.lines().collect();
    let from = lines
        .iter()
        .position(|l| l.trim_start().starts_with(anchor))
        .ok_or_else(|| format!("no `{anchor}` in the Rust source"))?;
    let at = from
        + lines[from..]
            .iter()
            .position(|l| l.trim_start().starts_with(item))
            .ok_or_else(|| format!("no `{item}` after `{anchor}` in the Rust source"))?;
    let mut docs = Vec::new();
    let mut k = at;
    while k > 0 {
        k -= 1;
        let l = lines[k].trim_start();
        if l.starts_with("#[") {
            continue;
        }
        match l.strip_prefix("///") {
            Some(rest) => docs.push(rest.strip_prefix(' ').unwrap_or(rest).to_string()),
            None => break,
        }
    }
    docs.reverse();
    let mut out = String::new();
    for d in docs {
        if d.is_empty() {
            out.push_str("#\n");
        } else {
            let _ = writeln!(out, "# {d}");
        }
    }
    Ok(out)
}

/// `line` as a `#` comment, with no trailing space when it is empty (lint-clean).
pub fn comment_line(line: &str) -> String {
    if line.is_empty() {
        "#".to_string()
    } else {
        format!("# {line}")
    }
}

/// The fixed header of every generated `.cyr` file.
pub fn header(title: &str, provenance: &str, notes: &[&str]) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "# tanmatra — {title}");
    out.push_str("#\n");
    let _ = writeln!(
        out,
        "# GENERATED by scripts/oracle from {provenance}. Do not edit."
    );
    out.push_str(
        "# Regenerate: cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml\n",
    );
    for n in notes {
        out.push_str("#\n");
        for line in n.lines() {
            let _ = writeln!(out, "# {line}");
        }
    }
    out.push('\n');
    out
}

/// A Cyrius array element type and how one value is written.
pub enum Elem {
    F64(Vec<f64>),
    U8(Vec<u8>),
    Bool(Vec<bool>),
    /// Names of constants (e.g. enum members), stored in a `u8` table.
    Names(Vec<String>),
}

/// `var NAME: T[N] = { .. };`, one element per line with a trailing label comment.
pub fn table(name: &str, doc: &str, elem: &Elem, labels: &[String]) -> Result<String, String> {
    let (ty, vals): (&str, Vec<String>) = match elem {
        Elem::F64(v) => (
            "f64",
            v.iter().map(|&x| f64_lit(x)).collect::<Result<_, _>>()?,
        ),
        Elem::U8(v) => ("u8", v.iter().map(u8::to_string).collect()),
        Elem::Bool(v) => ("bool", v.iter().map(bool::to_string).collect()),
        Elem::Names(v) => ("u8", v.clone()),
    };
    if vals.len() != labels.len() {
        return Err(format!(
            "{name}: {} values, {} labels",
            vals.len(),
            labels.len()
        ));
    }
    let mut out = String::new();
    for line in doc.lines() {
        let _ = writeln!(out, "{}", comment_line(line));
    }
    let _ = writeln!(out, "var {name}: {ty}[{}] = {{", vals.len());
    for (v, label) in vals.iter().zip(labels) {
        let _ = writeln!(out, "    {v},  # {label}");
    }
    out.push_str("};\n\n");
    Ok(out)
}
