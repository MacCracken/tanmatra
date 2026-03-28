# Architecture Overview

## Module Map

```
tanmatra (lib.rs)
  |-- constants.rs    CODATA 2022 fundamental constants
  |-- error.rs        TanmatraError enum
  |-- particle.rs     Standard Model: quarks, leptons, bosons, forces
  |-- nucleus.rs      Nuclear structure: binding energy, radius, magic numbers
  |-- decay.rs        Radioactive decay: modes, half-lives, chains
  |-- atomic.rs       Electron structure: configs, spectral lines, ionization
  |-- reaction.rs     Nuclear reactions: fusion, fission, Coulomb barrier
```

## Data Flow

```
Constants (CODATA 2022)
    |
    +-> Particle properties (masses, charges, spins)
    |
    +-> Nucleus (Z, A)
    |     |
    |     +-> Binding energy (Bethe-Weizsacker)
    |     +-> Nuclear radius
    |     +-> Mass defect
    |
    +-> Decay
    |     |
    |     +-> Decay constant (from half-life)
    |     +-> Activity (Bq)
    |     +-> Decay chains
    |
    +-> Atomic
    |     |
    |     +-> Electron configuration (Aufbau)
    |     +-> Spectral lines (Rydberg)
    |     +-> Ionization energies (NIST)
    |
    +-> Reactions
          |
          +-> Q-values
          +-> Coulomb barrier
```

## Consumers

- **prakash**: Optics and light simulation (spectral data)
- **kiran**: Game engine (physics simulation)
- **joshua**: Game manager (nuclear simulation scenarios)
- Any AGNOS component needing atomic/nuclear physics

## Design Decisions

- **Flat crate**: No workspace, single `src/` directory
- **`no_std` by default**: Core computations work without standard library
- **Real data only**: All constants from CODATA 2022, masses from PDG, ionization energies from NIST
- **Semi-empirical models**: Bethe-Weizsacker for binding energy (approximate but fast), Rydberg for spectral lines (exact for hydrogen-like)
- **Error handling**: All fallible operations return `Result<T, TanmatraError>`
