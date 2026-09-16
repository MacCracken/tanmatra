# ADR-005: Semi-Empirical Formulas with Corrections

**Date**: 2026-03-28 (revised 2026-09-16)
**Status**: Accepted (revised)

## Context

Nuclear binding energies can be computed from first principles (computationally expensive) or from semi-empirical formulas (fast, approximate). tanmatra needs to be fast enough for real-time use in kiran while maintaining physics accuracy.

## Decision

Use the Bethe-Weizsacker semi-empirical mass formula as the base, with coefficients least-squares fitted to the 2484 experimental AME2020 binding energies with A ≥ 16, and an optional shell correction in the published form of Myers & Swiatecki (Nucl. Phys. 81, 1 (1966)), whose amplitude C is fitted to the same data with the published c = 0.325.

Two methods: `binding_energy()` (RMS 3.31 MeV) and `binding_energy_shell_corrected()` (RMS 2.76 MeV; improves 1560 of 2484 nuclides; mean error for doubly magic nuclei −5.4 → −1.5 MeV).

The 2026-09-16 audit found that the previous "Strutinsky" correction (a sum of positive Gaussians with unsourced amplitudes) made results worse for 98.6% of nuclides, and that the textbook coefficients (fitted with a Z² Coulomb term but used with Z(Z−1)) gave RMS 10.2 MeV. Both were replaced.

## Consequences

**Positive**:
- `binding_energy()` is O(1) and suitable for hot loops
- Both formulas and their accuracy are traceable to AME2020
- Both methods available lets consumers choose speed vs accuracy

**Negative**:
- Semi-empirical formula is less accurate for light nuclei (A < 16)
- No deformation energy: deformed mid-shell nuclei retain errors of several MeV
