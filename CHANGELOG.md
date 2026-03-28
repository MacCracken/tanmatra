# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [0.1.0] - 2026-03-26

### Added

- `constants` module: CODATA 2022 fundamental constants (13 constants)
- `particle` module: Standard Model quarks, leptons, bosons, fundamental forces with real PDG masses
- `nucleus` module: Bethe-Weizsacker binding energy, mass defect, nuclear radius, magic numbers, 7 preset nuclei
- `decay` module: Decay modes, decay constant, activity, alpha/beta/gamma decay functions, 10 known isotopes with real half-lives, decay chain computation
- `atomic` module: Orbital types, quantum number validation, Aufbau electron configuration (with Cr/Cu exceptions), Rydberg spectral lines, Balmer series, NIST ionization energies for Z=1-36
- `reaction` module: Nuclear reactions with Q-values, D-T/D-D fusion, p-p chain, U-235 fission, Coulomb barrier
- `error` module: Comprehensive error types with thiserror
- `prelude` module: Convenient re-exports
- `no_std` support (default `std` feature)
- Integration test suite with physics validation
- Criterion benchmarks: binding_energy, spectral_line, electron_config, decay_chain
- CI/CD: GitHub Actions (check, test, MSRV, coverage, release)
