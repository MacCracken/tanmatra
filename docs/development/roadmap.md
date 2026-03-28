# Development Roadmap

## Completed

### v0.1.0 — Foundation (2026-03-26)

- [x] CODATA 2022 physical constants
- [x] Standard Model particles (quarks, leptons, bosons, forces)
- [x] Bethe-Weizsacker binding energy formula
- [x] Nuclear radius, magic numbers, preset nuclei
- [x] Radioactive decay (alpha, beta+, beta-, decay chains)
- [x] 17 known isotopes with real NNDC half-lives
- [x] Electron configurations (Aufbau + 5 exceptions)
- [x] Rydberg spectral lines
- [x] NIST ionization energies (Z=1-36)
- [x] Nuclear reactions with Q-values and Coulomb barriers
- [x] 7 preset reactions (DT, DD, pp, CNO, triple-alpha, U-235 fission)
- [x] Full serde support
- [x] `no_std` + `alloc` support
- [x] Criterion benchmarks
- [x] CI/CD pipelines

## Backlog

### v0.2.0 — Extended Nuclear Data

- [ ] Extended isotope database (100+ isotopes)
- [ ] Nuclear shell model (Woods-Saxon potential)
- [ ] Ionization energies for Z=37-118
- [ ] Electron affinity data
- [ ] Nuclear isomers and metastable states

### v0.3.0 — Quantum Mechanics

- [ ] Hydrogen atom wavefunctions (radial + angular)
- [ ] Transition probabilities (Einstein A/B coefficients)
- [ ] Selection rules for spectral transitions
- [ ] Zeeman and Stark effects
- [ ] Fine structure corrections

### v0.4.0 — Advanced Nuclear

- [ ] Liquid drop model improvements (shell corrections)
- [ ] Nuclear cross-sections (Breit-Wigner)
- [ ] Neutron moderation and thermalization
- [ ] Fission product yield distributions
- [ ] Bateman equations for decay chains

## Future

- Relativistic kinematics (Lorentz transformations, invariant mass)
- Scattering theory (Rutherford, Mott)
- Nuclear astrophysics (r-process, s-process)
- QED corrections to spectral lines
- Integration with prakash (optics) for photon interactions

## v1.0 Criteria

- 200+ tests
- 90%+ coverage
- 10+ benchmarks
- Complete Standard Model catalog
- Ionization energies for all elements
- 50+ known isotopes with validated half-lives
- All spectral series (Lyman, Balmer, Paschen, Brackett, Pfund)
- Peer-reviewed accuracy against NIST/NNDC data
