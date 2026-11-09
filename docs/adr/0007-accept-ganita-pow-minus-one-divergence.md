# 0007 — Accept ganita's pow(−1, |y| > 2^64) divergence until ganita fixes it

**Status**: Accepted (temporary; superseded when ganita ships the fix)
**Date**: 2026-11-08

## Context

The parity harness replays `tests/golden/libm.tsv`, the Rust libm 0.2.16 reference for every math
function tanmatra calls (port plan §2 step 8). One row diverges beyond any tolerance:

| Row | Rust libm | ganita 1.2.15 (`ganita_f64_pow`) |
|---|---|---|
| `libm::pow(-1.0, 1.7976931348623157e308)` | 1.0 | +inf |

Every finite |y| > 2^64 is an even integer, so the exact answer is 1. ganita's negative-base branch
calls its fdlibm core with |x| = 1, and that core lacks fdlibm's |x| == 1 guard. The defect is
filed upstream with a self-proving repro:
`ganita/docs/development/issues/2026-11-08-pow-minus-one-huge-even-exponent.md`. The repro exits
12, the number of wrong rows, on 1.2.15.

tanmatra itself is unaffected. Its three `pow` call sites all have positive bases:
- `remaining_fraction`: `pow(0.5, t / T½)`, `decay.rs:82`;
- `atomic_mass_amu`'s electron binding: `pow(Z, 2.39)` and `pow(Z, 5.35)` with Z ≥ 1,
  `nucleus.rs:357`;
- `radial_dipole_integral`: `pow(2/n, i)`, `atomic.rs:941`.

The row exists only because `libm.tsv` probes every function on the full special-value grid.

## Decision

The parity harness lists this one row, by function name and exact input bits, as an **accepted
divergence**. It is reported as `accepted (ADR 0007)` on every run and doesn't fail the module.
No other row is covered, and the class-C bound for `libm::pow` is unchanged at 2 ulp.

## Consequences

- **Positive:** M0's gate can close on a real bound for `pow`, and the defect stays visible in
  every harness run.
- **Negative:** until ganita ships the fix, a tanmatra caller that passed base −1 with |y| > 2^64
  would get ∞ or 0. No tanmatra function does.
- **Neutral:** when the cyrius pin moves to a stdlib with the fix, the row comes back ok on its
  own. At that point, remove the entry from `_par_accepted` in `tests/parity.tcyr` and set this
  ADR's status to *Superseded*.

## Alternatives considered

- **Widen pow's bound until the row passes.** Rejected. The port plan §4 forbids widening a class,
  and no relative bound can cover ∞ against 1.
- **Wrap pow in a tanmatra guard for |x| == 1.** Rejected for now. No tanmatra caller needs it,
  and it would mask the upstream fix. If a later module ever takes a negative base, this is the
  fallback.
- **Drop the row from `libm.tsv`.** Rejected. The corpus records the oracle as it is. Removing an
  inconvenient row is exactly the silent widening the policy exists to prevent.
