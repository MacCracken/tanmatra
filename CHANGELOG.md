# Changelog

All notable changes to tanmatra will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **nucleus**: Nuclear shell model (Mayer-Jensen, 32 levels through 184), ground-state spin-parity, shell occupation, shell closure functions, Strutinsky shell correction to binding energy
- **decay**: Extended isotope database from 17 to 114 isotopes (complete Th-232, U-235, U-238 chains, medical, industrial, cosmogenic, fission products, actinides), nuclear isomers (Ta-180m, Hf-178m2, Tc-99m, Am-242m, Pa-234m), isomeric transition decay mode, Bateman equations for sequential decay chain populations
- **particle**: Decay widths and lifetimes for W, Z, Higgs bosons and muon, tau leptons (PDG 2024), `is_stable()` method
- **atomic**: Ionization energies extended to Z=1-118 (NIST + theoretical), electron affinity data for Z=1-118, 17 new electron configuration exceptions (22 total), hydrogen atom radial wavefunctions (n=1-3), radial probability densities, Einstein A/B coefficients, electric dipole selection rules, fine-structure energy levels, named spectral series (Lyman, Balmer, Paschen, Brackett, Pfund), Zeeman effect (anomalous, Lande g-factor), Stark effect (linear, hydrogen), Lamb shift, vacuum polarization
- **reaction**: Cross-sections (geometric, Breit-Wigner resonance, thermal neutron database), neutron moderation (lethargy, collisions to thermalize, moderating ratio), fission product yield distributions (U-235, Pu-239), nucleosynthesis pathways (s-process, r-process)
- **relativity**: New module with four-momentum type, Lorentz factor, velocity addition, relativistic energy-momentum relations, de Broglie wavelength, invariant mass calculations
- **scattering**: New module with Rutherford differential/total cross-sections, Mott correction, distance of closest approach, Sommerfeld parameter
- **optics**: New feature-gated module (`optics` feature) integrating with prakash crate for spectral line to SPD conversion, Balmer series generator, wavelength-to-RGB color mapping
- **constants**: Bohr magneton, reduced Planck constant in MeV*s
- 11 criterion benchmarks (up from 4)
- 242 tests with 91.7% coverage

### Changed

- Binding energy asymmetry coefficient corrected from 93.15 to 23.285 MeV (standard convention)
- Pairing term exponent changed from A^(3/4) to A^(1/2) for consistency with coefficient set
- Coulomb constant corrected to 1.439_964_5 MeV*fm (derived from alpha*hbar*c)
- Top quark mass updated to 172_570 MeV (PDG 2024)
- Tau mass updated to 1_776.93 MeV (PDG 2024)
- W boson mass updated to 80_369 MeV (PDG 2024)
- 3 ionization energies corrected (H, Zn, Br) to match NIST
- 5 isotope half-lives updated to NUBASE2020 (I-131, Cs-137, Rn-222, U-234, Am-241)
- `decay_chain()` performance improved 8x by hoisting isotope database allocation
- `#[inline]` added to all hot-path computation functions

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
