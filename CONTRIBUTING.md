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
- `#[inline]` on hot-path computation functions
- Serde (`Serialize + Deserialize`) on all public types
- Zero `unwrap`/`panic` in library code
- `no_std` compatible -- use `alloc` not `std` collections
- All physics constants from CODATA 2022
- All ionization energies from NIST

## License

By contributing, you agree that your contributions will be licensed under GPL-3.0-only.
