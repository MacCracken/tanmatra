# ADR-005: Semi-Empirical Formulas with Corrections

**Date**: 2026-03-28
**Status**: Accepted

## Context

Nuclear binding energies can be computed from first principles (computationally expensive) or from semi-empirical formulas (fast, approximate). tanmatra needs to be fast enough for real-time use in kiran while maintaining physics accuracy.

## Decision

Use the Bethe-Weizsacker semi-empirical mass formula as the base, with optional Strutinsky shell corrections for improved accuracy near magic numbers. The shell model (Mayer-Jensen) provides the correction data.

Two methods: `binding_energy()` (fast, ~1-2% error) and `binding_energy_shell_corrected()` (slower, better near magic numbers).

## Consequences

**Positive**:
- `binding_energy()` is O(1) and suitable for hot loops
- Shell correction improves accuracy for doubly-magic nuclei
- Both methods available lets consumers choose speed vs accuracy

**Negative**:
- Semi-empirical formula is less accurate for light nuclei (A < 20)
- Strutinsky correction is parameterized, not derived from first principles
