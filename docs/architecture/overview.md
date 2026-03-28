# Architecture Overview

## Module Map

```
tanmatra/
  src/
    lib.rs          — crate root, prelude, feature flags
    error.rs        — TanmatraError (6 variants, thiserror)
    constants.rs    — CODATA 2022 physical constants
    particle.rs     — Standard Model: quarks, leptons, bosons, forces
    nucleus.rs      — Nuclear structure, Bethe-Weizsacker binding energy
    decay.rs        — Radioactive decay, half-lives, decay chains
    atomic.rs       — Electron configuration, spectral lines, ionization
    reaction.rs     — Nuclear reactions, Q-values, Coulomb barriers
  benches/
    benchmarks.rs   — Criterion benchmarks (4 suites)
  examples/
    basic.rs        — Demonstrates all modules
```

## Data Flow

```
constants.rs ──> All modules (physical constants)
      |
particle.rs     (standalone — Standard Model catalog)
      |
nucleus.rs ───> decay.rs ───> decay chains
      |              |
      +──────> reaction.rs ──> Q-values, barriers
      |
atomic.rs       (standalone — electron structure)
```

## Dependencies

- **libm**: `no_std` math functions (cbrt, pow, log)
- **serde**: Serialization for all public types
- **thiserror**: Error derive macro
- **tracing** (optional): Structured logging

## Design Decisions

1. **Flat crate**: Single-level modules, no nested crate hierarchy
2. **Real data only**: All constants from CODATA 2022, masses from PDG 2024, half-lives from NNDC
3. **`no_std` first**: Core functionality works without `std`
4. **Semi-empirical formulas**: Bethe-Weizsacker for binding energy (accurate to ~1% for A > 20)
5. **Known exceptions**: Electron configuration handles Cr, Cu, Mo, Ag, Au anomalies
