# ADR-004: Feature-Gated Optics Integration

**Date**: 2026-03-28
**Status**: Accepted

## Context

tanmatra computes spectral line wavelengths. prakash (AGNOS optics crate) can visualize these as spectral power distributions. We needed to decide how to connect them without creating a hard dependency.

## Decision

The `optics` module is gated behind an `optics` feature flag that pulls in `prakash` as an optional dependency. tanmatra defines a `SpectralLine` struct and provides conversion functions to prakash's `Spd` type.

## Consequences

**Positive**:
- No dependency overhead for consumers who don't need optics
- Clean separation: tanmatra computes, prakash visualizes
- `SpectralLine` is useful even without prakash (wavelength + intensity)

**Negative**:
- Two feature flags to remember (`optics` implies `std`)
- API surface increases when feature is enabled
