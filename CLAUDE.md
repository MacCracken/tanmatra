# tanmatra — Claude Code Instructions

## Project Identity

**tanmatra** (Sanskrit: subtle element) — Atomic and subatomic physics for AGNOS

- **Type**: Flat library crate
- **License**: GPL-3.0
- **MSRV**: 1.89
- **Version**: SemVer 0.1.0

## Consumers

prakash (optics/light), kiran (game engine), joshua (game manager), and any AGNOS component needing physics computations.

## Development Process

### P(-1): Scaffold Hardening (before any new features)

0. Read roadmap, CHANGELOG, and open issues
1. Test + benchmark sweep of existing code
2. Cleanliness check: `cargo fmt --check`, `cargo clippy --all-features --all-targets -- -D warnings`, `cargo audit`, `cargo deny check`, `RUSTDOCFLAGS="-D warnings" cargo doc --all-features --no-deps`
3. Get baseline benchmarks
4. Internal deep review
5. External research — CODATA, NIST, PDG, nuclear data tables
6. Cleanliness check — must be clean after review
7. Additional tests/benchmarks from findings
8. Post-review benchmarks
9. Repeat if heavy

### Work Loop (continuous)

1. Work phase
2. Cleanliness check
3. Test + benchmark additions
4. Run benchmarks
5. Internal review
6. Cleanliness check
7. Deeper tests/benchmarks
8. Benchmarks again
9. If review heavy -> return to step 5
10. Documentation — CHANGELOG, roadmap, docs
11. Version check
12. Return to step 1

### Key Principles

- Never skip benchmarks
- `#[non_exhaustive]` on ALL public enums
- `#[must_use]` on all pure functions
- Every type must be Serialize + Deserialize (serde)
- Feature-gate optional modules
- Zero unwrap/panic in library code
- All types must have serde roundtrip tests
- ALL physics values must use real data (CODATA 2022, NIST, PDG, NNDC)
- No fake data, no stubs, no placeholder values

## DO NOT

- **Do not commit or push** — the user handles all git operations
- **NEVER use `gh` CLI** — use `curl` to GitHub API only
- Do not add unnecessary dependencies
- Do not break backward compatibility without a major version bump
- Do not skip benchmarks before claiming performance improvements
- Do not use approximate/made-up physics values — always cite sources
