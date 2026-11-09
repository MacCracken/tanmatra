//! `oracle` — writes everything the Cyrius port derives from the frozen Rust tanmatra 1.3.0
//! in `rust-old/`:
//! - the golden corpus, `tests/golden/<module>.tsv` (port plan §2 step 5);
//! - the generated data modules, `src/constants.cyr` and `src/<module>_data.cyr`
//!   (port plan §2 step 6).
//!
//! ```sh
//! cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml            # write
//! cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check # verify
//! ```
//!
//! `--check` regenerates every output in memory and compares it byte for byte with the
//! committed file. It exits 1 if any output differs or is missing, or a `.tsv` in
//! `tests/golden/` has no emitter.
//!
//! Modules are added at the start of their milestone, with the golden grid from the port
//! plan §3. M1 so far: constants, particle, relativity; plus `libm`, the math-function
//! reference for the M0 step 8 baselines.

mod cyr;
mod grids;
mod m_constants;
mod m_libm;
mod m_particle;
mod m_relativity;
mod t_constants;
mod t_particle;
mod tsv;

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

/// Every golden emitter, in file-name order.
fn golden_modules() -> Vec<tsv::Module> {
    vec![
        m_constants::emit(),
        m_libm::emit(),
        m_particle::emit(),
        m_relativity::emit(),
    ]
}

/// One generated file: its path from the repo root, its contents, and a summary.
struct Output {
    rel: String,
    text: String,
    summary: String,
}

fn read(path: &Path) -> Result<String, String> {
    fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))
}

/// Every generated Cyrius source, in file-name order. Each is built from the oracle's API.
/// The Rust source is read for doc comments only.
fn cyrius_outputs(root: &Path, prov: &str) -> Result<Vec<Output>, String> {
    let rust = |m: &str| read(&root.join(format!("rust-old/src/{m}.rs")));
    let mut out = Vec::new();
    let mut push = |rel: &str, text: String| {
        let lines = text.lines().count();
        out.push(Output {
            rel: rel.to_string(),
            text,
            summary: format!("{lines} lines"),
        });
    };
    push(
        "src/constants.cyr",
        t_constants::emit(&rust("constants")?, prov)?,
    );
    push(
        "src/particle_data.cyr",
        t_particle::emit(&rust("particle")?, prov)?,
    );
    Ok(out)
}

/// `version = "x"` from the first `[package]` line that has it.
fn toml_version(text: &str) -> Option<&str> {
    let line = text.lines().find(|l| l.starts_with("version = \""))?;
    line.strip_prefix("version = \"")?.strip_suffix('"')
}

/// The locked version of `crate_name` in a Cargo.lock.
fn lock_version<'a>(lock: &'a str, crate_name: &str) -> Option<&'a str> {
    let needle = format!("name = \"{crate_name}\"\n");
    let at = lock.find(&needle)? + needle.len();
    let line = lock[at..].lines().next()?;
    line.strip_prefix("version = \"")?.strip_suffix('"')
}

fn provenance(root: &Path) -> Result<String, String> {
    let cargo = read(&root.join("rust-old/Cargo.toml"))?;
    let lock = read(&root.join("scripts/oracle/Cargo.lock"))?;
    let tm = toml_version(&cargo).ok_or("rust-old/Cargo.toml: no version")?;
    let lm = lock_version(&lock, "libm").ok_or("scripts/oracle/Cargo.lock: no libm")?;
    Ok(format!("tanmatra {tm} (rust-old/), libm {lm}"))
}

/// Write `text` to `root/rel`, or in check mode compare it. Returns false on a mismatch.
fn emit(root: &Path, check: bool, out: &Output) -> Result<bool, String> {
    let path: PathBuf = root.join(&out.rel);
    if !check {
        fs::write(&path, &out.text).map_err(|e| format!("{}: {e}", path.display()))?;
        println!("wrote    {} ({})", out.rel, out.summary);
        return Ok(true);
    }
    match fs::read_to_string(&path) {
        Ok(old) if old == out.text => {
            println!("ok       {} ({})", out.rel, out.summary);
            Ok(true)
        }
        Ok(_) => {
            println!("DIFFERS  {}", out.rel);
            Ok(false)
        }
        Err(_) => {
            println!("MISSING  {}", out.rel);
            Ok(false)
        }
    }
}

fn run() -> Result<bool, String> {
    let check = std::env::args().skip(1).any(|a| a == "--check");
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let root = root
        .canonicalize()
        .map_err(|e| format!("{}: {e}", root.display()))?;
    let prov = provenance(&root)?;
    let dir = root.join("tests/golden");
    if !check {
        fs::create_dir_all(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    }

    let mut outputs = Vec::new();
    let mut emitted = BTreeSet::new();
    for module in golden_modules() {
        let name = format!("{}.tsv", module.name);
        emitted.insert(name.clone());
        outputs.push(Output {
            rel: format!("tests/golden/{name}"),
            text: module.render(&prov),
            summary: format!("{} rows", module.len()),
        });
    }
    // The harness reads this list, so a module added here can't go untested: without a
    // dispatcher it is reported as pending, row by row.
    let mut index = String::from(
        "# tanmatra golden corpus: one module per line; tests/parity.tcyr replays each\n",
    );
    for name in &emitted {
        index.push_str(name.trim_end_matches(".tsv"));
        index.push('\n');
    }
    outputs.push(Output {
        rel: "tests/golden/modules.txt".to_string(),
        summary: format!("{} modules", emitted.len()),
        text: index,
    });
    outputs.extend(cyrius_outputs(&root, &prov)?);

    let mut ok = true;
    for out in &outputs {
        ok &= emit(&root, check, out)?;
    }

    // A .tsv with no emitter is stale: the harness would test rows nobody regenerates.
    if let Ok(entries) = fs::read_dir(&dir) {
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().into_owned();
            if name.ends_with(".tsv") && !emitted.contains(&name) {
                println!("STALE    tests/golden/{name} (no emitter)");
                ok = false;
            }
        }
    }
    Ok(ok)
}

fn main() -> ExitCode {
    match run() {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::from(1),
        Err(e) => {
            eprintln!("oracle: {e}");
            ExitCode::from(2)
        }
    }
}
