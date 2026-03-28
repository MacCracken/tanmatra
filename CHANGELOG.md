# Changelog

All notable changes to tanmatra will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-03-26

### Added

- **constants**: CODATA 2022 physical constants (electron/proton/neutron masses, fine-structure constant, Rydberg constant, Bohr radius, and more)
- **particle**: Standard Model implementation with quarks (6 flavors, PDG 2024 masses), leptons (6 types), gauge bosons (photon, gluon, W, Z, Higgs), and four fundamental forces with relative strengths
- **nucleus**: Nuclear structure with Bethe-Weizsacker semi-empirical mass formula (binding energy, B/A, mass defect, nuclear radius), magic number identification, and presets (H-1, He-4, C-12, Fe-56, U-235, U-238)
- **decay**: Radioactive decay calculations (decay constant, remaining fraction, activity in Bq), alpha/beta-minus/beta-plus decay transformations, 17 known isotopes with real NNDC half-lives (H-3 through Am-241), and decay chain generation
- **atomic**: Quantum number validation, electron configurations via Aufbau/Madelung with exceptions (Cr, Cu, Mo, Ag, Au), Rydberg spectral line wavelengths, NIST ionization energies for Z=1-36, noble gas core notation
- **reaction**: Nuclear reaction framework with Q-value computation, Coulomb barrier estimation, and 7 preset reactions (DT fusion, DD fusion x2, pp chain, U-235 fission, CNO cycle, triple-alpha)
- **error**: TanmatraError with 6 variants (InvalidAtomicNumber, InvalidMassNumber, InvalidQuantumNumbers, InvalidHalfLife, DecayNotPossible, InvalidReaction)
- Full serde support on all public types
- `no_std` + `alloc` support
- Criterion benchmarks (binding_energy, spectral_line, electron_config, decay_chain)
- CI/CD pipelines (check, test, coverage, MSRV, release)
