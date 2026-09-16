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
Accuracies below were measured against reference data (see
`docs/audit/2026-09-16-math-audit.md`).

| Computation | Accuracy | Limitation |
|-------------|----------|------------|
| Binding energy (Bethe-Weizsacker, AME2020 fit) | RMS 3.31 MeV over 2484 nuclides (A ≥ 16) | Liquid-drop model; no deformation |
| Binding energy (+ Myers–Swiatecki shell term) | RMS 2.76 MeV | Spherical shell term only |
| Spectral lines, `spectral_line_nm` | Exact Rydberg formula for infinite nuclear mass | Real lines are longer by 1 + m_e/M (5.4e-4 for H); use `spectral_line_vacuum_nm` |
| Fine-structure energy (first order) | Agrees with the Dirac energy to O((Zα)⁴) | No QED, no reduced mass |
| Einstein A (hydrogen-like) | Exact non-relativistic dipole rates; with μ/m_e matches NIST to < 2e-4 | No relativistic corrections |
| Lamb shift (one loop, Zα expansion) | ≈0.2% for H, ≈3% for He⁺ | Not valid for high Z; n ≤ 4 |
| Rutherford / Mott scattering | Exact for point charges (first Born for Mott) | No screening, recoil or finite size unless requested |
| Pair production (Bethe–Heitler, Maximon) | Analytic unscreened Born result | No screening (≈30% high for Pb at 1 GeV), no Coulomb correction |
| Bateman chains | Relative accuracy ~1e-12 for arbitrary (including equal) decay constants | Matrix-exponential cost O(n³ log(λt)) |
| `AtomicInstant` arithmetic | Exact to 1 ns | Conversions to/from f64 seconds are limited by f64 resolution |
| Julian dates (f64) | ≈40 µs resolution | Use `AtomicInstant` for sub-millisecond work |
