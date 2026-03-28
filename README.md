# tanmatra

**tanmatra** (Sanskrit: तन्मात्र — subtle element) — Atomic and subatomic physics for [AGNOS](https://github.com/MacCracken/agnosticos).

Standard Model particles, nuclear structure (Bethe-Weizsacker), radioactive decay with real half-lives, spectral lines (Rydberg), electron configurations (Aufbau + exceptions), nuclear reactions, and ionization energies.

## Features

- **Standard Model**: Quarks, leptons, bosons with PDG 2024 masses; four fundamental forces with coupling strengths
- **Nuclear structure**: Binding energy (semi-empirical mass formula), nuclear radii, magic numbers
- **Radioactive decay**: Decay modes, decay constants, activity, 17 known isotopes with NNDC half-lives, decay chains
- **Atomic physics**: Quantum number validation, electron configurations with Cr/Cu/Mo/Ag/Au exceptions, Rydberg spectral lines
- **Nuclear reactions**: Q-values, Coulomb barriers, preset fusion (DT, DD, pp, CNO, triple-alpha) and fission (U-235)
- **Ionization energies**: NIST values for Z=1-36
- **`no_std` compatible**: Works with `alloc` only, no standard library required
- **Serde support**: All types serialize/deserialize

## Quick Start

```rust
use tanmatra::prelude::*;

// Iron-56 binding energy per nucleon (~8.8 MeV)
let fe56 = Nucleus::iron_56();
println!("B/A = {:.2} MeV", fe56.binding_energy_per_nucleon());

// H-alpha spectral line (~656.3 nm)
let h_alpha = spectral_line_nm(1, 2, 3).unwrap();
println!("H-alpha = {:.1} nm", h_alpha);

// Electron configuration of iron: [Ar] 4s2 3d6
let config = electron_configuration(26).unwrap();
println!("{}", format_configuration_short(&config, 26));
```

## Data Sources

| Source | Used For |
|--------|----------|
| [CODATA 2022](https://physics.nist.gov/cuu/Constants/) | Fundamental constants |
| [PDG 2024](https://pdg.lbl.gov/) | Particle masses |
| [NNDC/NUBASE](https://www.nndc.bnl.gov/) | Nuclear half-lives |
| [NIST ASD](https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html) | Ionization energies |

## Feature Flags

| Feature | Default | Description |
|---------|---------|-------------|
| `std` | Yes | Standard library support |
| `logging` | No | Structured tracing |
| `full` | No | All optional features |

## License

GPL-3.0-only. See [LICENSE](LICENSE).
