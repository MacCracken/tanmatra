# ADR-003: `no_std` First

**Date**: 2026-03-26
**Status**: Accepted

## Context

tanmatra may be used in embedded or WASM contexts (via kiran game engine). The library needs to work without the standard library.

## Decision

The crate is `no_std` by default with `alloc` for heap allocation. The `std` feature (enabled by default) adds standard library support. All math uses `libm` instead of `std::f64` methods.

## Consequences

**Positive**:
- Works in embedded, WASM, and kernel contexts
- Forces explicit allocation awareness
- `libm` provides consistent cross-platform math

**Negative**:
- Cannot use `std::collections::HashMap` without `std` feature
- `libm` functions may be slightly slower than hardware-optimized `std` equivalents
