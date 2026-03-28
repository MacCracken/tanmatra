# Security Policy

## Scope

tanmatra is a pure physics computation library. It performs no I/O, no network access, and contains no `unsafe` code. All computations are deterministic.

## Attack Surface

| Area | Risk | Mitigation |
|------|------|------------|
| Binding energy formula | Division by zero for A=0 | Validated in constructor, returns error |
| Decay calculations | Division by zero for zero half-life | Half-life validated > 0 |
| Quantum numbers | Invalid combinations | Validated in constructor |
| Serde deserialization | Crafted JSON with extreme values | Enum validation via serde derive; parameters validated on construction |
| Nuclear radius | Negative mass number | Validated in constructor |

## Reporting Vulnerabilities

Report security issues to the repository maintainer via GitHub Security Advisories. Do not file public issues for security vulnerabilities.

## Dependencies

| Dependency | Purpose | Risk |
|---|---|---|
| `serde` | Serialization | Widely audited, no unsafe in derive |
| `thiserror` | Error derive | Proc macro only, no runtime code |
| `libm` | `no_std` math | Pure Rust, no unsafe |
| `tracing` (optional) | Structured logging | No I/O; subscriber is caller's responsibility |
