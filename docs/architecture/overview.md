# Architecture Overview

## Module Map

```
tanmatra/
  src/
    lib.rs          -- crate root, prelude, feature flags, Send+Sync assertions
    error.rs        -- TanmatraError (6 variants, thiserror)
    constants.rs    -- CODATA 2022 physical constants (17 constants)
    particle.rs     -- Standard Model: quarks, leptons, bosons, forces (PDG 2024)
    nucleus.rs      -- Bethe-Weizsacker, shell model, Strutinsky correction
    decay.rs        -- Radioactive decay, 114 isotopes, Bateman equations
    atomic.rs       -- Electron config, spectral series, Zeeman/Stark, QED, wavefunctions
    reaction.rs     -- Q-values, cross-sections, fission yields, moderation, nucleosynthesis
    relativity.rs   -- Four-momentum, Lorentz transformations, de Broglie
    scattering.rs   -- Rutherford/Mott scattering cross-sections
    optics.rs       -- prakash integration (feature-gated: "optics")
  benches/
    benchmarks.rs   -- Criterion benchmarks (11 suites)
  examples/
    basic.rs        -- Core module demonstration
    nuclear.rs      -- Nuclear physics: binding energy, shell model, decay chains
    spectral.rs     -- Spectral series and fine structure
    relativity.rs   -- Relativistic kinematics
  tests/
    integration.rs  -- Cross-module integration tests (20 tests)
```

## Data Flow

```
constants.rs ──────────> All modules (physical constants)
      |
particle.rs              Standalone: Standard Model catalog, decay widths
      |
nucleus.rs ────────────> decay.rs ──────> decay chains, Bateman equations
      |  (binding energy)    |  (isotopes)
      |  shell model         |
      +────────────────> reaction.rs ──> Q-values, cross-sections, nucleosynthesis
      |                      |
      |                 fission yields, neutron moderation
      |
atomic.rs ──────────────> optics.rs (feature-gated)
  (spectral lines)          (prakash SPD conversion)
  (wavefunctions)
  (Zeeman/Stark/QED)

relativity.rs            Standalone: four-vectors, Lorentz kinematics
scattering.rs            Uses: constants (fine-structure, hbar*c)
```

## Consumers

| Crate | Uses |
|-------|------|
| [prakash](https://crates.io/crates/prakash) | Spectral lines for optical simulation |
| [kiran](https://crates.io/crates/kiran) | Physics computations for game engine |
| Any AGNOS component | Nuclear, atomic, or particle physics |

## Dependencies

| Crate | Purpose | Optional |
|-------|---------|----------|
| `libm` | `no_std` math (sqrt, cbrt, exp, log, sin, cos) | No |
| `serde` | Serialization for all public types | No |
| `thiserror` | Error derive macro | No |
| `tracing` | Structured logging | Yes (`logging`) |
| `prakash` | Optics/spectral integration | Yes (`optics`) |

## Design Decisions

See [docs/decisions/](decisions/) for full ADRs.

1. **Flat crate** (ADR-001): Single-level modules, no nested hierarchy
2. **Real data only** (ADR-002): All values from CODATA, PDG, NIST, NNDC
3. **`no_std` first** (ADR-003): Core functionality works without `std`
4. **Feature-gated optics** (ADR-004): prakash integration behind `optics` feature
5. **Semi-empirical + corrections** (ADR-005): Bethe-Weizsacker + Strutinsky shell correction

## Public Types

| Module | Types | Traits |
|--------|-------|--------|
| `particle` | `Quark`, `Lepton`, `Boson`, `FundamentalForce` | Serialize, Deserialize, Copy, Hash |
| `nucleus` | `Nucleus`, `ShellLevel` | Serialize, Deserialize, Copy, Hash |
| `decay` | `DecayMode`, `Isotope` | Serialize, Deserialize, Clone |
| `atomic` | `OrbitalType`, `QuantumNumbers`, `OrbitalFilling`, `TransitionType` | Serialize, Deserialize |
| `reaction` | `NuclearReaction`, `ThermalCrossSection`, `FissionYield`, `NucleosynthesisProcess`, `NucleosynthesisStep`, `NucleosynthesisPathway` | Serialize, Deserialize |
| `relativity` | `FourMomentum` | Serialize, Deserialize, Copy, Add |
| `error` | `TanmatraError` | Serialize, Deserialize, thiserror::Error |
