//! Strict extraction of the oracle's PRIVATE literal constants and tables.
//!
//! Some values the port needs have no public accessor: the SEMF and Myers–Swiatecki
//! coefficients, AME2020's u→keV factor, and the empirical proton/neutron level orders. They
//! must not be re-typed (port plan §2 step 6), so they are read mechanically from
//! `rust-old/src`, under these rules:
//! - Only a `const NAME: TYPE = <literal>;` item is read. A literal is a decimal number, or for
//!   a table, a list of integer tuples. Anything else is an error, never a guess.
//! - A table's element count must equal its declared length.
//! - A float is parsed by `str::parse::<f64>`. That is correctly rounded, as rustc's literal
//!   conversion is, so the bits are the oracle's.
//! - Every extracted value is then PROVEN against the oracle's public behaviour by its emitter:
//!   it recomputes something the oracle computes from the value, and requires identical bits.
//!   An extraction that slipped would fail there, before anything is written.

/// The text of `const NAME: ... = <body>;` (the body only), or an error.
fn const_body<'a>(src: &'a str, name: &str) -> Result<&'a str, String> {
    let key = format!("const {name}:");
    let at = src
        .find(&key)
        .ok_or_else(|| format!("private const {name} not found"))?;
    if src[at + key.len()..].contains(&key) {
        return Err(format!("private const {name} is defined twice"));
    }
    let eq = src[at..]
        .find('=')
        .ok_or_else(|| format!("{name}: no `=`"))?
        + at;
    let end = src[eq..]
        .find(';')
        .ok_or_else(|| format!("{name}: no `;`"))?
        + eq;
    Ok(src[eq + 1..end].trim())
}

/// A private `const NAME: f64 = <decimal literal>;`.
pub fn f64_const(src: &str, name: &str) -> Result<f64, String> {
    let body = const_body(src, name)?;
    let lit = body.replace('_', "");
    let ok = lit
        .chars()
        .all(|c| c.is_ascii_digit() || matches!(c, '.' | '-' | 'e' | 'E' | '+'));
    if !ok || lit.is_empty() {
        return Err(format!("{name}: `{body}` is not a plain decimal literal"));
    }
    lit.parse::<f64>()
        .map_err(|e| format!("{name}: `{body}`: {e}"))
}

/// A private `const NAME: [u32; LEN] = [a, b, ...];`.
pub fn u32_array(src: &str, name: &str, len: usize) -> Result<Vec<u32>, String> {
    let body = const_body(src, name)?;
    let inner = body
        .strip_prefix('[')
        .and_then(|b| b.strip_suffix(']'))
        .ok_or_else(|| format!("{name}: not an array literal"))?;
    let vals = inner
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<u32>().map_err(|e| format!("{name}: `{s}`: {e}")))
        .collect::<Result<Vec<_>, _>>()?;
    if vals.len() != len {
        return Err(format!("{name}: {} elements, declared {len}", vals.len()));
    }
    Ok(vals)
}

/// A private `const NAME: [(u32, u32, u32); LEN] = [(a, b, c), ...];`. Comments are skipped.
pub fn u32_triples(src: &str, name: &str, len: usize) -> Result<Vec<(u32, u32, u32)>, String> {
    let body = const_body(src, name)?;
    let text: String = body
        .lines()
        .map(|l| l.split("//").next().unwrap_or(""))
        .collect::<Vec<_>>()
        .join(" ");
    let mut out = Vec::new();
    let mut rest = text.as_str();
    while let Some(open) = rest.find('(') {
        let close = rest[open..]
            .find(')')
            .ok_or_else(|| format!("{name}: unclosed tuple"))?
            + open;
        let parts = rest[open + 1..close]
            .split(',')
            .map(|s| {
                s.trim()
                    .parse::<u32>()
                    .map_err(|e| format!("{name}: `{s}`: {e}"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        if parts.len() != 3 {
            return Err(format!("{name}: a tuple of {} elements", parts.len()));
        }
        out.push((parts[0], parts[1], parts[2]));
        rest = &rest[close + 1..];
    }
    if out.len() != len {
        return Err(format!("{name}: {} tuples, declared {len}", out.len()));
    }
    Ok(out)
}

/// A private `const NAME: &[(f64, ..., f64)] = &[(a, b, ...), ...];` slice of `arity`-tuples of
/// decimal literals. A slice declares no length, so the count is proven by the emitter's
/// recomputation, never assumed here.
pub fn f64_tuples(src: &str, name: &str, arity: usize) -> Result<Vec<Vec<f64>>, String> {
    let body = const_body(src, name)?;
    let text: String = body
        .lines()
        .map(|l| l.split("//").next().unwrap_or(""))
        .collect::<Vec<_>>()
        .join(" ");
    let inner = text
        .trim()
        .strip_prefix("&[")
        .and_then(|b| b.trim_end().strip_suffix(']'))
        .ok_or_else(|| format!("{name}: not a slice literal"))?;
    let mut out = Vec::new();
    let mut rest = inner;
    while let Some(open) = rest.find('(') {
        let close = rest[open..]
            .find(')')
            .ok_or_else(|| format!("{name}: unclosed tuple"))?
            + open;
        let parts = rest[open + 1..close]
            .split(',')
            .map(|s| {
                let lit = s.trim().replace('_', "");
                let plain = !lit.is_empty()
                    && lit
                        .chars()
                        .all(|c| c.is_ascii_digit() || matches!(c, '.' | '-' | 'e' | 'E' | '+'));
                if !plain {
                    return Err(format!(
                        "{name}: `{}` is not a plain decimal literal",
                        s.trim()
                    ));
                }
                lit.parse::<f64>()
                    .map_err(|e| format!("{name}: `{lit}`: {e}"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        if parts.len() != arity {
            return Err(format!(
                "{name}: a tuple of {} elements, expected {arity}",
                parts.len()
            ));
        }
        out.push(parts);
        rest = &rest[close + 1..];
    }
    if out.is_empty() {
        return Err(format!("{name}: no tuples"));
    }
    Ok(out)
}
