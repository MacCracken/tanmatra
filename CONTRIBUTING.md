# Contributing to tanmatra

Thank you for your interest in contributing to tanmatra.

## Development Workflow

1. Fork and clone the repository
2. Create a feature branch from `main`
3. Make your changes
4. Run the cleanliness check (see below)
5. Open a pull request

## Prerequisites

- The Cyrius toolchain at the version pinned in `cyrius.cyml` (`[package].cyrius`)
- Rust stable (MSRV 1.89), only to run the frozen oracle in `rust-old/` and the generators in `scripts/oracle/`

tanmatra is being ported from Rust to Cyrius. The plan of record is
[`docs/development/tanmatra-port-plan.md`](docs/development/tanmatra-port-plan.md).
`rust-old/` is the parity oracle: never edit it.

## Cleanliness Check

Every change must pass:

```bash
cyrius fmt --check src/*.cyr tests/*.tcyr
cyrius lint src/*.cyr
cyrius build src/main.cyr build/tanmatra
cyrius test            # must also be warning-free
cyrius test --bounds
cyrius bench
cyrius deps --locked
cyrius audit
```

## Code Conventions

- Public enums have explicit, append-only numeric values; never renumber them
- `#must_use` on all pure functions
- No aborts in library code: no `+?` / `-?` / `*?`, `assert_fatal` or `result_unwrap` in `src/`
- f64 hygiene: bind call results to `: f64` locals; no hex in a `const`; `-x`, never `(0 - x)`, on an f64
- All physics values must use CODATA 2022 / NIST / PDG / NNDC data. Tables are generated, never typed by hand
- Feature-gate optional modules (`[features]` + `[lib.<profile>]`)

## Porting a Module

1. Read the module's milestone in the port plan (§3) and its golden grid
2. Write `src/<module>.cyr` and append it to `[lib] modules` in dependency order
3. Add its golden rows (`scripts/oracle/`) and port its reference-value tests to `tests/*.tcyr`
4. Add a bench to `tests/*.bcyr` and record it against the Rust number
5. Check surface parity: every in-scope Rust `pub fn` has a Cyrius counterpart or a written reason
6. Update the CHANGELOG, `docs/development/state.md` and the port plan's tracking table

## License

By contributing, you agree that your contributions will be licensed under GPL-3.0-only.
