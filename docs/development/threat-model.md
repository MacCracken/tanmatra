# Threat Model

## Trust Boundaries

tanmatra is a pure computation library. It does not perform I/O, network access, or file operations. All inputs are numeric values passed by the caller.

```
Untrusted ──> [Input validation] ──> [Pure computation] ──> [Result]
              (Z, A, n, l, energy)    (no I/O, no unsafe)    (f64 or Error)
```

## Attack Surface

| Area | Risk | Mitigation |
|------|------|------------|
| Integer overflow in Z, A, n | Incorrect physics results | All checked at API boundary; `u32` prevents negative values |
| Floating-point edge cases | NaN/Inf propagation | Division-by-zero guarded; `is_finite()` checks on inputs |
| Large allocations | Memory exhaustion | `known_isotopes()` is bounded (114 entries); `decay_chain()` has `max_steps` limit |
| Serde deserialization | Invalid state construction | `Nucleus` fields are private; validation via `new()` constructor |
| Supply chain (dependencies) | Malicious code | 4 runtime deps, all audited; `cargo deny` + `cargo audit` in CI |
| Unsafe code | Memory safety | `#![forbid(unsafe_code)]` at crate root |

## Panic Sites

**Production code: ZERO.** The crate uses `#![forbid(unsafe_code)]` and has no `unwrap()`, `expect()`, `panic!()`, `unreachable!()`, or `todo!()` in non-test code.

`unwrap_or_else` is used in `known_isotopes()` for infallible `Nucleus::new()` calls where inputs are compile-time constants. The fallback values are never reached.

## Supply Chain

| Dependency | Version | Purpose | Risk |
|------------|---------|---------|------|
| `libm` | 0.2 | `no_std` math | Low: pure Rust, no deps, widely used |
| `serde` | 1 | Serialization | Low: ubiquitous, audited |
| `thiserror` | 2 | Error derives | Low: proc-macro only, no runtime code |
| `tracing` | 0.1 | Logging (optional) | Low: widely audited |
| `prakash` | 1.1 | Optics (optional) | Low: AGNOS ecosystem, same author |

## Numeric Precision

All computations use `f64` (IEEE 754 double precision, ~15-16 significant digits).

| Computation | Typical Error | Source of Error |
|-------------|--------------|-----------------|
| Binding energy (Bethe-Weizsacker) | ~1-2% for A > 20 | Semi-empirical formula limitations |
| Binding energy (shell-corrected) | ~0.5-1% near magic numbers | Parameterized correction |
| Spectral lines (Rydberg) | < 0.01% for hydrogen | Formula is exact for hydrogen-like |
| Fine-structure correction | ~0.001% | First-order perturbation theory |
| Rutherford scattering | Exact (classical limit) | Formula is exact for point Coulomb |
| Bateman equations | < 0.01% | Numerical precision of exp() |
