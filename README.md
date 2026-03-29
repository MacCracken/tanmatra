# tanmatra

**tanmatra** (Sanskrit: तन्मात्र — subtle element) — Atomic and subatomic physics for [AGNOS](https://github.com/MacCracken/agnosticos).

[![CI](https://github.com/MacCracken/tanmatra/actions/workflows/ci.yml/badge.svg)](https://github.com/MacCracken/tanmatra/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/tanmatra.svg)](https://crates.io/crates/tanmatra)
[![docs.rs](https://docs.rs/tanmatra/badge.svg)](https://docs.rs/tanmatra)
[![license](https://img.shields.io/crates/l/tanmatra.svg)](LICENSE)

## Modules

| Module | Description |
|--------|-------------|
| `constants` | CODATA 2022 fundamental physical constants |
| `particle` | Standard Model: quarks, leptons, bosons, forces (PDG 2024) |
| `nucleus` | Bethe-Weizsacker binding energy, shell model, magic numbers |
| `decay` | Radioactive decay, 114 isotopes, Bateman equations, decay chains |
| `atomic` | Electron configs, spectral series, ionization/affinity, Zeeman/Stark, QED |
| `reaction` | Q-values, cross-sections, fission yields, neutron moderation, nucleosynthesis |
| `relativity` | Four-momentum, Lorentz factor, velocity addition, de Broglie |
| `scattering` | Rutherford/Mott cross-sections, Sommerfeld parameter |
| `optics` | prakash integration for spectral line visualization (feature-gated) |
| `error` | `TanmatraError` with 6 variants |

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

// Relativistic proton at 500 MeV/c
let p = FourMomentum::from_mass_and_momentum(PROTON_MASS_MEV, 500.0);
println!("gamma = {:.3}", p.gamma());
```

## Building

```bash
make check    # fmt + clippy + test + audit
make bench    # criterion benchmarks (11 suites)
make doc      # rustdoc with -D warnings
make coverage # tarpaulin coverage report
```

**MSRV**: Rust 1.89

## Data Sources

| Source | Used For |
|--------|----------|
| [CODATA 2022](https://physics.nist.gov/cuu/Constants/) | Fundamental constants |
| [PDG 2024](https://pdg.lbl.gov/) | Particle masses, decay widths |
| [NNDC/NUBASE 2020](https://www.nndc.bnl.gov/) | Nuclear half-lives (114 isotopes) |
| [NIST ASD](https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html) | Ionization energies, electron affinities |
| [ENDF/B-VIII.0](https://www.nndc.bnl.gov/endf/) | Cross-sections, fission yields |

## Feature Flags

| Feature | Default | Description |
|---------|---------|-------------|
| `std` | Yes | Standard library support |
| `logging` | No | Structured tracing via `tracing` crate |
| `optics` | No | Integration with [prakash](https://crates.io/crates/prakash) for spectral visualization |
| `full` | No | Enables all optional features |

## License

GPL-3.0-only. See [LICENSE](LICENSE).
