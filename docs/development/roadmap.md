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

### v0.2.0 — Extended Nuclear Data (2026-03-28)

- [x] Extended isotope database (114 isotopes with NNDC half-lives)
- [x] Complete Th-232, U-235, U-238 natural decay chains
- [x] Nuclear shell model (Mayer-Jensen, 32 levels, spin-parity)
- [x] Strutinsky shell correction to binding energy
- [x] Nuclear isomers (Ta-180m, Hf-178m2, Tc-99m, Am-242m)
- [x] Ionization energies for Z=1-118 (NIST + theoretical)
- [x] Electron affinity data for Z=1-118
- [x] Electron configuration exceptions (22 total, all periods)
- [x] Particle decay widths and lifetimes (PDG 2024)
- [x] PDG 2024 mass updates (top, tau, W)
- [x] Corrected CODATA constants (Coulomb constant, asymmetry coefficient)

### v0.3.0 — Quantum Mechanics (2026-03-28)

- [x] Hydrogen atom wavefunctions (radial R_nl for n=1-3)
- [x] Radial probability densities
- [x] Transition probabilities (Einstein A and B coefficients)
- [x] Electric dipole selection rules (Δl, Δm_l)
- [x] Fine-structure corrections (Dirac energy levels)
- [x] Zeeman effect (anomalous, Lande g-factor)
- [x] Stark effect (linear, hydrogen)
- [x] Named spectral series (Lyman, Balmer, Paschen, Brackett, Pfund)

### v0.4.0 — Advanced Nuclear (2026-03-28)

- [x] Liquid drop model improvements (Strutinsky shell corrections)
- [x] Nuclear cross-sections (geometric, Breit-Wigner, thermal neutron)
- [x] Thermal neutron cross-section database (8 key isotopes)
- [x] Neutron moderation and thermalization (lethargy, collisions)
- [x] Fission product yield distributions (U-235, Pu-239)
- [x] Bateman equations for sequential decay chain populations
- [x] QED corrections (Lamb shift, vacuum polarization)
- [x] Nucleosynthesis pathways (s-process, r-process)

### v0.5.0 — Relativity and Scattering (2026-03-28)

- [x] Relativistic kinematics (four-momentum, Lorentz factor, velocity addition)
- [x] Energy-momentum relations, de Broglie wavelength
- [x] Invariant mass calculations
- [x] Rutherford scattering (differential, total, distance of closest approach)
- [x] Mott scattering correction
- [x] Sommerfeld parameter
- [x] Integration with prakash (optics) — spectral line to SPD conversion

## v1.0 Criteria — Status

- [x] 200+ tests (242 tests)
- [x] 90%+ coverage (91.7%)
- [x] 10+ benchmarks (11 benchmarks)
- [x] Complete Standard Model catalog (PDG 2024)
- [x] Ionization energies for all elements (Z=1-118)
- [x] 50+ known isotopes (114 isotopes)
- [x] All spectral series (Lyman, Balmer, Paschen, Brackett, Pfund)
- [x] Peer-reviewed accuracy against NIST/NNDC data
