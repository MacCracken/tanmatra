//! `golden` — writes the Cyrius port's golden corpus, `tests/golden/<module>.tsv`, from the
//! frozen Rust tanmatra 1.3.0 in `rust-old/` (port plan §2 step 5).
//!
//! ```sh
//! cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml            # write
//! cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check # verify
//! ```
//!
//! `--check` regenerates every module in memory and compares it byte for byte with the
//! committed file. It exits 1 if any module differs, is missing, or a stale `.tsv` has no
//! emitter.
//!
//! Modules are added at the start of their milestone, with the golden grid from the port
//! plan §3. M1 so far: constants, particle, relativity; plus `libm`, the math-function
//! reference for the M0 step 8 baselines.

mod grids;
mod m_constants;
mod m_libm;
mod m_particle;
mod m_relativity;
mod tsv;

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

/// Every emitter, in file-name order.
fn modules() -> Vec<tsv::Module> {
    vec![
        m_constants::emit(),
        m_libm::emit(),
        m_particle::emit(),
        m_relativity::emit(),
    ]
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
    let read = |p: PathBuf| fs::read_to_string(&p).map_err(|e| format!("{}: {e}", p.display()));
    let cargo = read(root.join("rust-old/Cargo.toml"))?;
    let lock = read(root.join("scripts/oracle/Cargo.lock"))?;
    let tm = toml_version(&cargo).ok_or("rust-old/Cargo.toml: no version")?;
    let lm = lock_version(&lock, "libm").ok_or("scripts/oracle/Cargo.lock: no libm")?;
    Ok(format!("tanmatra {tm} (rust-old/), libm {lm}"))
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

    let mut ok = true;
    let mut emitted = BTreeSet::new();
    for module in modules() {
        let name = format!("{}.tsv", module.name);
        let path = dir.join(&name);
        let text = module.render(&prov);
        emitted.insert(name.clone());
        if check {
            match fs::read_to_string(&path) {
                Ok(old) if old == text => println!("ok       {name} ({} rows)", module.len()),
                Ok(_) => {
                    println!("DIFFERS  {name}");
                    ok = false;
                }
                Err(_) => {
                    println!("MISSING  {name}");
                    ok = false;
                }
            }
        } else {
            fs::write(&path, &text).map_err(|e| format!("{}: {e}", path.display()))?;
            println!("wrote    {name} ({} rows)", module.len());
        }
    }

    // A .tsv with no emitter is stale: the harness would test rows nobody regenerates.
    if let Ok(entries) = fs::read_dir(&dir) {
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().into_owned();
            if name.ends_with(".tsv") && !emitted.contains(&name) {
                println!("STALE    {name} (no emitter)");
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
            eprintln!("golden: {e}");
            ExitCode::from(2)
        }
    }
}
