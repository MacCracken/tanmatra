# ADR-001: Flat Module Layout

**Date**: 2026-03-26
**Status**: Accepted

## Context

tanmatra covers multiple physics domains (particles, nuclei, atoms, reactions, relativity, scattering). We needed to decide between a nested crate workspace or a flat single-crate layout.

## Decision

Use a flat library crate with one module per physics domain. All modules live in `src/` at the top level.

## Consequences

**Positive**:
- Simple dependency graph (no inter-crate versioning)
- Easy cross-module access via `crate::` paths
- Single `Cargo.toml` to maintain
- Consumers add one dependency

**Negative**:
- Compile times grow with crate size
- Feature-gating is the only way to exclude modules
