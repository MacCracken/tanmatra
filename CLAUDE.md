# tanmatra — Claude Code Instructions

> **Core rule**: this file is **preferences, process, and procedures** —
> durable rules that change rarely. Volatile state (current version,
> module line counts, port progress, test counts, consumers) lives in
> [`docs/development/state.md`](docs/development/state.md).
> Do not inline state here.

## Project Identity

**tanmatra** (Sanskrit: subtle element) — Atomic and subatomic physics for AGNOS

- **Type**: Flat library (Cyrius); a Rust → Cyrius port in progress
- **License**: GPL-3.0-only
- **Language**: Cyrius. The toolchain is pinned in `cyrius.cyml` `[package].cyrius`
- **Version**: `VERSION` at the project root is the source of truth. Don't inline the number here
- **Rust oracle**: tanmatra 1.3.0, frozen at `rust-old/`. It is the parity reference until the port completes
- **Plan of record**: [`docs/development/tanmatra-port-plan.md`](docs/development/tanmatra-port-plan.md)
- **Standards**: [First-Party Standards](https://github.com/MacCracken/agnosticos/blob/main/docs/development/applications/first-party-standards.md) · [First-Party Documentation](https://github.com/MacCracken/agnosticos/blob/main/docs/development/applications/first-party-documentation.md)

## Goal

tanmatra owns atomic, nuclear and particle physics for AGNOS: Standard Model data,
nuclear structure and decay, atomic spectra, reactions, scattering and atomic
timekeeping. Every value comes from an evaluated primary source.

## Consumers

prakash (optics/light), kiran (game engine), joshua (game manager), and any AGNOS component needing physics computations.

## Quick Start

```sh
cyrius deps                                   # resolve dependencies (vendors lib/)
cyrius build src/main.cyr build/tanmatra      # compile
cyrius test                                   # run tests/*.tcyr (parity + reference tests)
cyrius bench                                  # run tests/*.bcyr
cargo test --manifest-path rust-old/Cargo.toml --all-features   # the Rust oracle, read-only
```

## Development Process

### P(-1): Scaffold Hardening (before any new features)

0. Read the roadmap, the port plan, the CHANGELOG and the open issues
1. Test + benchmark sweep of the existing code
2. Cleanliness check (below)
3. Get baseline benchmarks
4. Internal deep review
5. External research — CODATA, NIST, PDG, nuclear data tables
6. Cleanliness check — must be clean after the review
7. Additional tests/benchmarks from the findings
8. Post-review benchmarks
9. Repeat if heavy

### Cleanliness check

```sh
cyrius fmt --check
cyrius lint
cyrius vet
cyrius build src/main.cyr build/tanmatra
cyrius test                 # must be warning-free, not just passing
cyrius test --bounds
cyrius bench
cyrius distlib --all && cyrius distlib --check    # once [lib] has modules
cyrius deps --locked
cyrius audit
```

`cyrius test` must be **warning-free**. The f64 type diagnostics come only from the
compiles that `test` and `build` drive, never from `cyrius lint`. Each one is a bug:
- *integer arithmetic with an f64 right operand*
- *f64 arithmetic with a non-f64 right operand*
- *an integer stored into an f64 slot*
- *comparison mixes f64 and integer*

### Work Loop (continuous)

1. Work phase — one module at a time, in the order the port plan gives
2. Cleanliness check
3. Test + benchmark additions: golden parity rows, reference-value tests, benches
4. Run benchmarks, recorded against the Rust numbers
5. Internal review, including a surface-parity diff against `rust-old/src/`
6. Cleanliness check
7. Deeper tests/benchmarks
8. Benchmarks again
9. If the review was heavy, return to step 5
10. Documentation — CHANGELOG, `state.md`, roadmap, the port plan's tracking table
11. Version check
12. Return to step 1

## Key Principles

- **Cross-check against `rust-old/`.** The port's correctness bar is "matches what Rust
  did", within the tolerance class the port plan (§4) assigns. Diverge only with an ADR
  in `docs/adr/`.
- **Correctness over cleverness.** If the Cyrius behavior diverges silently from Rust, the bugs win.
- Test after every change, not after the feature is "done".
- ONE change at a time — never bundle unrelated changes.
- Never skip benchmarks.
- **Enums.** Public enums have explicit, append-only numeric values that are never
  renumbered. Property functions `match` with an `_ =>` arm.
- **`#must_use` on all pure functions**, plus `#pure` where it applies.
- **JSON (bayan / `#derive(Serialize)`) only for the exchange types** listed in the port plan (M7).
- **Feature-gate optional modules** with `[features]` plus a `[lib.<profile>]` bundle.
- **No aborts in library code.** `src/` never uses `+?`, `-?` or `*?` (they exit the
  process), `assert_fatal` or `result_unwrap`. Errors return a value-form `Result` or an
  `err_out` code.
- **f64 hygiene.**
  - Bind every call result to a `: f64` local before using it in an operator.
  - Never put hex in a `const`: it is an integer.
  - Never write `(0 - x)` on an f64: use `-x`.
  - Use the canonical `ganita_f64_*` names, not the deprecated short aliases.
- **Tables and constants are generated, never typed by hand.** The generators live in
  `scripts/oracle/`.
- Build with `cyrius build`, not raw `cat file | cc5` — the manifest auto-resolves deps.
- Source files only need project includes — stdlib auto-resolves from `cyrius.cyml`.
- `var buf[N]` = N **bytes**, not N entries.
- ALL physics values must use real data (CODATA 2022, NIST, PDG, NNDC).
- No fake data, no stubs, no placeholder values.

## DO NOT

- **Do not commit or push** — the user handles all git operations
- **NEVER use `gh` CLI** — use `curl` to the GitHub API only
- **Do not modify `rust-old/`** — it is the parity oracle. Tools that drive it (golden and table generators) live in `scripts/oracle/` and depend on it by path
- Do not modify `lib/` files (vendored stdlib / dep symlinks)
- Do not hardcode toolchain versions in CI YAML — `cyrius = "X.Y.Z"` in `cyrius.cyml` is the source of truth
- Do not widen a parity tolerance class to make a row pass
- Do not add unnecessary dependencies
- Do not break backward compatibility without a major version bump
- Do not skip benchmarks before claiming performance improvements
- Do not skip tests before claiming changes work
- Do not use approximate/made-up physics values — always cite sources

## Documentation

- [`docs/development/tanmatra-port-plan.md`](docs/development/tanmatra-port-plan.md) — the port plan of record
- [`docs/development/state.md`](docs/development/state.md) — live state
- [`docs/development/roadmap.md`](docs/development/roadmap.md) — milestones and backlog
- [`docs/adr/`](docs/adr/) — Architecture Decision Records (*why X over Y?*)
- [`docs/architecture/`](docs/architecture/) — non-obvious constraints
- [`docs/guides/`](docs/guides/) — task-oriented how-tos
- [`docs/audit/`](docs/audit/) — audits (the 2026-09-16 math audit is the correctness baseline)
- `rust-old/docs/` — the Rust 1.3.0 testing guide and architecture overview, kept for reference
