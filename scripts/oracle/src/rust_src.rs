//! Strict extraction of the oracle's PRIVATE literal constants and tables.
//!
//! Some values the port needs have no public accessor: the SEMF and Myers–Swiatecki
//! coefficients, AME2020's u→keV factor, the empirical proton/neutron level orders, the
//! Aufbau filling order and the Bethe logarithms. They must not be re-typed (port plan §2
//! step 6), so they are read mechanically from `rust-old/src`, under these rules:
//! - Only a `const NAME: TYPE = <literal>;` item is read. A literal is a decimal number, or for
//!   a table, a list of integer tuples (or `(integer, Enum::Variant)` pairs). The one other
//!   shape is a pure lookup function, `Some(match (a, b) { (1, 0) => <decimal>, .. _ => return
//!   None, })`, every line of which must be such an arm. Anything else is an error, never a
//!   guess.
//! - Comments are read for labels only (`z_comments`: each table entry's citation).
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
    // The terminating `;` is the first one outside a `//` comment (citations contain `;`).
    let mut end = None;
    let mut in_comment = false;
    let bytes = src.as_bytes();
    for k in eq..bytes.len() {
        match bytes[k] {
            b'\n' => in_comment = false,
            b'/' if !in_comment && bytes.get(k + 1) == Some(&b'/') => in_comment = true,
            b';' if !in_comment => {
                end = Some(k);
                break;
            }
            _ => {}
        }
    }
    let end = end.ok_or_else(|| format!("{name}: no `;`"))?;
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

/// A private `const NAME: [(u32, ENUM); LEN] = [(a, ENUM::Variant), ...];`, as (integer,
/// variant name) pairs. Comments are skipped.
pub fn u32_variant_pairs(
    src: &str,
    name: &str,
    enum_name: &str,
    len: usize,
) -> Result<Vec<(u32, String)>, String> {
    let body = const_body(src, name)?;
    let text: String = body
        .lines()
        .map(|l| l.split("//").next().unwrap_or(""))
        .collect::<Vec<_>>()
        .join(" ");
    let prefix = format!("{enum_name}::");
    let mut out = Vec::new();
    let mut rest = text.as_str();
    while let Some(open) = rest.find('(') {
        let close = rest[open..]
            .find(')')
            .ok_or_else(|| format!("{name}: unclosed tuple"))?
            + open;
        let parts: Vec<&str> = rest[open + 1..close].split(',').map(str::trim).collect();
        if parts.len() != 2 {
            return Err(format!("{name}: a tuple of {} elements", parts.len()));
        }
        let n = parts[0]
            .parse::<u32>()
            .map_err(|e| format!("{name}: `{}`: {e}", parts[0]))?;
        let v = parts[1]
            .strip_prefix(&prefix)
            .filter(|v| !v.is_empty() && v.chars().all(|c| c.is_ascii_alphanumeric()))
            .ok_or_else(|| format!("{name}: `{}` is not a {enum_name} variant", parts[1]))?;
        out.push((n, v.to_string()));
        rest = &rest[close + 1..];
    }
    if out.len() != len {
        return Err(format!("{name}: {} tuples, declared {len}", out.len()));
    }
    Ok(out)
}

/// One lookup arm: the `(a, b)` pattern and its value.
pub type Arm = ((u32, u32), f64);

/// The arms `(a, b) => <decimal literal>,` of a private lookup function
/// `fn NAME(..) -> Option<f64> { Some(match (..) { .. _ => return None, }) }`. Every line of
/// the match must be such an arm, except its `_ => return None,` fallback; anything else is
/// an error. Like a const table, the values must be proven by the emitter's recomputation.
pub fn match_f64_arms(src: &str, fn_name: &str) -> Result<Vec<Arm>, String> {
    let key = format!("fn {fn_name}(");
    let at = src
        .find(&key)
        .ok_or_else(|| format!("private fn {fn_name} not found"))?;
    if src[at + key.len()..].contains(&key) {
        return Err(format!("private fn {fn_name} is defined twice"));
    }
    let body = &src[at..];
    let open = body
        .find("Some(match")
        .ok_or_else(|| format!("{fn_name}: not a `Some(match ..)` lookup"))?;
    let mut out = Vec::new();
    for line in body[open..].lines().skip(1) {
        let l = line.trim();
        if l == "_ => return None," {
            return Ok(out);
        }
        let (lhs, rhs) = l
            .split_once("=>")
            .ok_or_else(|| format!("{fn_name}: `{l}` is not a match arm"))?;
        let pair = lhs
            .trim()
            .strip_prefix('(')
            .and_then(|p| p.strip_suffix(')'))
            .ok_or_else(|| format!("{fn_name}: `{l}`: the pattern is not a pair"))?;
        let (a, b) = pair
            .split_once(',')
            .ok_or_else(|| format!("{fn_name}: `{l}`: the pattern is not a pair"))?;
        let a = a
            .trim()
            .parse::<u32>()
            .map_err(|e| format!("{fn_name}: `{l}`: {e}"))?;
        let b = b
            .trim()
            .parse::<u32>()
            .map_err(|e| format!("{fn_name}: `{l}`: {e}"))?;
        let lit = rhs
            .trim()
            .strip_suffix(',')
            .ok_or_else(|| format!("{fn_name}: `{l}`: no trailing comma"))?
            .replace('_', "");
        let plain = !lit.is_empty()
            && lit
                .chars()
                .all(|c| c.is_ascii_digit() || matches!(c, '.' | '-' | 'e' | 'E' | '+'));
        if !plain {
            return Err(format!("{fn_name}: `{l}`: not a plain decimal literal"));
        }
        let v = lit
            .parse::<f64>()
            .map_err(|e| format!("{fn_name}: `{lit}`: {e}"))?;
        out.push(((a, b), v));
    }
    Err(format!("{fn_name}: no `_ => return None,` arm"))
}

/// The trailing `//` comment of each element line of a private `const NAME: [..; LEN]` table,
/// in order, for labels only (values always come from the oracle's API). Each must name its
/// element as `(Z=k)`, k = 1, 2, ..; a table out of step with Z is an error.
pub fn z_comments(src: &str, name: &str, len: usize) -> Result<Vec<String>, String> {
    let body = const_body(src, name)?;
    let mut out = Vec::new();
    for line in body.lines() {
        if let Some((_, c)) = line.split_once("//") {
            let c = c.trim().to_string();
            let want = format!("(Z={})", out.len() + 1);
            if !c.contains(&want) {
                return Err(format!("{name}: comment `{c}` does not name {want}"));
            }
            out.push(c);
        }
    }
    if out.len() != len {
        return Err(format!(
            "{name}: {} commented entries, declared {len}",
            out.len()
        ));
    }
    Ok(out)
}
