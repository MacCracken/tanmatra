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
| `nucleus` | AME2020-fitted binding energy + Myers–Swiatecki shell term, shell model, AME2020 masses, radii, moments, superallowed Ft |
| `decay` | Radioactive decay, 114 isotopes (NUBASE2020), stable Bateman solver, decay chains |
| `atomic` | Electron configs, spectral lines, ionization/affinity, exact hydrogenic wavefunctions and A-values, Zeeman/Stark, Lamb shift, Dirac |
| `reaction` | AME2020 Q-values, ENDF/B-VIII.0 thermal cross sections, resonance integrals, fission yields, moderation, nucleosynthesis |
| `relativity` | Four-momentum, Lorentz factor, velocity addition, de Broglie |
| `scattering` | Rutherford/Mott/Born cross-sections, form factors, Klein–Nishina, Bethe–Heitler pair production |
| `timekeeping` | CIPM 2025 frequency standards, TAI/UTC/TT/GPS/TCG/TCB/TDB, leap seconds, relativistic clock corrections |
| `optics` | prakash integration for spectral line visualization (feature-gated) |
| `error` | `TanmatraError` with 6 variants |

## Quick Start

```rust
use tanmatra::prelude::*;

// Iron-56 binding energy per nucleon (~8.8 MeV)
let fe56 = Nucleus::iron_56();
println!("B/A = {:.2} MeV", fe56.binding_energy_per_nucleon());

// H-alpha vacuum wavelength with the proton reduced mass (656.470 nm)
let h_alpha = spectral_line_vacuum_nm(1, 1, 2, 3).unwrap();
println!("H-alpha = {:.3} nm", h_alpha);

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
make bench    # criterion benchmarks (19 suites)
make doc      # rustdoc with -D warnings
make coverage # tarpaulin coverage report
```

**MSRV**: Rust 1.89

## Data Sources

| Source | Used For |
|--------|----------|
| [CODATA 2022](https://physics.nist.gov/cuu/Constants/) | Fundamental constants |
| [PDG 2024](https://pdg.lbl.gov/) | Particle masses, decay widths |
| [AME2020 / NUBASE2020](https://www-nds.iaea.org/amdc/) | Atomic masses; half-lives, decay modes, J^π (114 isotopes) |
| [NIST ASD](https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html) | Ionization energies (Z ≤ 108), hydrogen A-values |
| Andersen, Haugen & Hotop (1999) + later measurements; Smits et al. (2023) | Electron affinities; superheavy ionization energies |
| [ENDF/B-VIII.0](https://www.nndc.bnl.gov/endf/) | Thermal cross sections, resonance integrals, fission yields |
| Angeli & Marinova (2013); Stone INDC(NDS)-0794/0833 | Charge radii; magnetic dipole and quadrupole moments |
| Hardy & Towner (2020) | Superallowed ft and Ft values |
| CIPM 2025, IERS Conventions 2010, IERS Bulletin C | Frequency standards, time scales, leap seconds |
| Storey & Hummer (1995) | Case B Balmer intensities |

## Feature Flags

| Feature | Default | Description |
|---------|---------|-------------|
| `std` | Yes | Standard library support |
| `logging` | No | Structured tracing via `tracing` crate |
| `optics` | No | Integration with [prakash](https://crates.io/crates/prakash) for spectral visualization |
| `full` | No | Enables all optional features |

## License

GPL-3.0-only. See [LICENSE](LICENSE).
