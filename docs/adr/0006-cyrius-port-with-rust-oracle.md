# 0006 — Port to Cyrius, with Rust 1.3.0 frozen as the parity oracle

**Status**: Accepted
**Date**: 2026-11-08

## Context

AGNOS libraries are moving from Rust to Cyrius, and tanmatra's in-Cyrius consumer
(prakash) has already ported. The 2026-09-16 math audit found that the Rust crate was
wrong in 36 formulas and about 100 table values. Every finding was repaired in Rust
1.3.0, pinned by reference-value tests against primary sources. A translation needs an
oracle it can check against value by value, and only a repaired crate can be one.

The Rust crate was `no_std`-first and used `libm` for portable math (ADR 0003). Neither
mechanism exists in Cyrius.

## Decision

tanmatra 2.0.0 is a Cyrius library. The Rust 1.3.0 crate is frozen in `rust-old/` as the
parity oracle, and nothing in it is edited.
- **Initialization.** The port was started with `cyrius port` on Cyrius 6.7.18.
- **Scope.** The milestones, the scope (physics first; bridges and the soorat integration
  deferred until their consumers port) and the parity policy are those of
  `docs/development/tanmatra-port-plan.md`.
- **Parity classes.** Each golden value has a class: A exact, B arithmetic, C
  transcendental with a bound measured per target, or D iterative.
- **Divergences.** A divergence from the oracle beyond its class is a defect, unless an ADR
  accepts it.

## Consequences

- **Positive**
  - Correctness is measured, not asserted: every in-scope public function is diffed
    against the repaired Rust crate, bit pattern by bit pattern.
  - The audit's reference-value tests travel with the code.
- **Negative**
  - Two toolchains are needed until the port completes: cargo for the oracle and its
    generators, cyrius for the library.
  - `rust-old/` stays in the tree until every deferred bridge has ported or been dropped.
  - Cyrius exp/ln/pow/trig are within 1 ulp but not bit-identical to `libm`, and exp/ln
    differ between x86_64 and aarch64, so class C needs measured, per-target tolerances.
- **Neutral**
  - The generators that drive the oracle live in `scripts/oracle/` and depend on
    `rust-old/` by path, so the oracle itself is never modified.
  - Cargo features become `[features]` plus `[lib.<profile>]` bundles, and serde becomes
    JSON for the exchange types only.

## Alternatives considered

- **Port and fix together, with no frozen oracle.** Rejected. Without a fixed reference,
  a translation error and a deliberate fix look the same, and the audit showed that
  hand-carried numbers fail at a high rate.
- **Keep tanmatra in Rust.** Rejected. Its Cyrius consumers would need a foreign-function
  bridge for every physics call, against the AGNOS direction.
- **Translate from the pre-audit crate.** Rejected. Parity tests would bless the 36 wrong
  formulas.
