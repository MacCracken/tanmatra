# Contributing to tanmatra

Thank you for your interest in contributing to tanmatra.

## Development Workflow

1. Fork and clone the repository
2. Create a feature branch from `main`
3. Make your changes
4. Run the cleanliness check (see below)
5. Open a pull request

## Prerequisites

- Rust stable (MSRV 1.89)
- Components: `rustfmt`, `clippy`
- Optional: `cargo-audit`, `cargo-deny`

## Cleanliness Check

Every change must pass:

```bash
cargo fmt --check
cargo clippy --all-features --all-targets -- -D warnings
cargo test --all-features
cargo test --no-default-features
RUSTDOCFLAGS="-D warnings" cargo doc --all-features --no-deps
cargo audit
cargo deny check
```

## Code Conventions

- `#[non_exhaustive]` on all public enums
- `#[must_use]` on all pure functions
- Serde (`Serialize + Deserialize`) on all public types
- Zero `unwrap`/`panic` in library code
- `no_std` compatible -- use `alloc` not `std` collections
- All physics values must use CODATA 2022 / NIST data
- Feature-gate optional modules

## Adding a New Module

1. Create `src/my_module.rs` following the pattern in existing modules
2. Register in `lib.rs`: module declaration, prelude export, Send+Sync assertion
3. Add integration tests: physics accuracy, edge cases, serde roundtrip
4. Add a criterion benchmark
5. Document data sources (CODATA, NIST, etc.)

## License

By contributing, you agree that your contributions will be licensed under GPL-3.0-only.
