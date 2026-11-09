# Architecture Decision Records

Decisions about tanmatra — what we chose, the context, and the consequences we accept. Use these when a future reader would reasonably ask *"why did we do it this way?"*

## Conventions

- **Filename**: `NNNN-kebab-case-title.md`, zero-padded to four digits. Never renumber.
- **One decision per ADR.** If a decision supersedes a prior one, add a new ADR and set the old one's status to `Superseded by NNNN`.
- **Status lifecycle**: `Proposed` → `Accepted` → (optionally) `Superseded` or `Deprecated`.
- Use [`template.md`](template.md) as the starting point.

## ADR vs. architecture note vs. guide

| Kind | Lives in | Answers |
|---|---|---|
| ADR | `docs/adr/` | *Why did we choose X over Y?* |
| Architecture note | `docs/architecture/` | *What non-obvious constraint is true about the code?* |
| Guide | `docs/guides/` | *How do I do X?* |

## Index

| ADR | Title | Status |
|---|---|---|
| [0001](0001-flat-module-layout.md) | Flat module layout | Accepted |
| [0002](0002-real-data-only.md) | Real data only | Accepted |
| [0003](0003-no-std-first.md) | `no_std` first | Superseded by 0006 |
| [0004](0004-feature-gated-optics.md) | Feature-gated optics integration | Accepted |
| [0005](0005-semi-empirical-plus-corrections.md) | Semi-empirical formulas with corrections | Accepted (revised) |
| [0006](0006-cyrius-port-with-rust-oracle.md) | Port to Cyrius, with Rust 1.3.0 frozen as the parity oracle | Accepted |
| [0007](0007-accept-ganita-pow-minus-one-divergence.md) | Accept ganita's pow(−1, \|y\| > 2^64) divergence until ganita fixes it | Accepted (temporary) |

ADRs 0001–0005 were written for the Rust crate (as `docs/decisions/001`–`005`, renamed to the four-digit convention in 2026-11). Their decisions still hold for the port, except 0003.
