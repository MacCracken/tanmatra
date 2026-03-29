# ADR-002: Real Data Only

**Date**: 2026-03-26
**Status**: Accepted

## Context

Physics libraries can use approximate/textbook values or real evaluated data. For a library consumed by production simulation code (kiran, prakash), accuracy matters.

## Decision

All physics values must come from authoritative evaluated sources:
- Constants: CODATA 2022
- Particle masses: PDG 2024
- Half-lives: NNDC / NUBASE 2020
- Ionization energies: NIST ASD
- Cross-sections: ENDF/B-VIII.0

No fake data, no stubs, no placeholder values.

## Consequences

**Positive**:
- Downstream consumers can trust computed values
- Errors are traceable to specific source publications
- Results are reproducible against reference data

**Negative**:
- Data entry is tedious and error-prone (mitigated by verification tests)
- Superheavy element data (Z > 104) relies on theoretical predictions
