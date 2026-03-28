# Security Policy

## Scope

tanmatra is a pure physics computation library. It performs no I/O, no network access, and contains no `unsafe` code. All computations are deterministic.

## Attack Surface

| Area | Risk | Mitigation |
|------|------|------------|
| Binding energy formula | NaN/Infinity on edge inputs | Validated Z and A ranges |
| Spectral line computation | Division by zero on equal quantum numbers | Returns error for n_upper <= n_lower |
| Decay chain computation | Infinite loop on stable isotope | Step limit parameter required |
| Electron configuration | Stack overflow on extreme Z | Z limited to valid range (1-118) |
| Serde deserialization | Crafted JSON with extreme values | Enum validation via serde derive; parameters validated on use |

## Reporting Vulnerabilities

Report security issues to the repository maintainer via GitHub Security Advisories. Do not file public issues for security vulnerabilities.

## Dependencies

| Dependency | Purpose | Risk |
|---|---|---|
| `serde` | Serialization | Widely audited, no unsafe in derive |
| `thiserror` | Error derive | Proc macro only, no runtime code |
| `libm` | `no_std` math | Pure Rust, no unsafe |
| `tracing` (optional) | Structured logging | No I/O; subscriber is caller's responsibility |
