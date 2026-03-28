# Development Roadmap

## Completed

### v0.1.0 (2026-03-26)

- [x] CODATA 2022 fundamental constants
- [x] Standard Model particles (quarks, leptons, bosons, forces)
- [x] Nuclear binding energy (Bethe-Weizsacker)
- [x] Nuclear radius and magic numbers
- [x] Radioactive decay (alpha, beta, gamma, EC, fission)
- [x] 10 known isotopes with real half-lives
- [x] Decay chain computation
- [x] Electron configuration with Aufbau rule + exceptions
- [x] Spectral lines (Rydberg formula)
- [x] Balmer series
- [x] NIST ionization energies (Z=1-36)
- [x] Nuclear reactions (fusion, fission, Q-values)
- [x] Coulomb barrier
- [x] Comprehensive test suite
- [x] Criterion benchmarks
- [x] CI/CD pipelines

## Backlog

- [ ] Shell model energy levels
- [ ] Liquid drop model refinements (Strutinsky correction)
- [ ] More decay chains (real U-238 chain with mixed modes)
- [ ] Cross-section calculations
- [ ] Ionization energies for Z=37-118
- [ ] Electron configuration exceptions beyond Cr/Cu
- [ ] Relativistic corrections to spectral lines
- [ ] Zeeman/Stark effect
- [ ] Nuclear isomers
- [ ] Nucleosynthesis pathways (r-process, s-process)
- [ ] Particle decay widths and lifetimes

## Future

- [ ] Integration with prakash for spectral emission/absorption
- [ ] Monte Carlo decay simulation
- [ ] Nuclear data tables (AME2020)
- [ ] ENDF cross-section data integration
- [ ] Feynman diagram computation helpers

## v1.0 Criteria

- All Standard Model particles with PDG 2024 properties
- Binding energies within 1% of experimental for Z > 20
- Complete ionization energies (Z=1-118)
- At least 50 known isotopes
- Full decay chain for U-238 and Th-232 series
- >90% test coverage
- Published on crates.io
