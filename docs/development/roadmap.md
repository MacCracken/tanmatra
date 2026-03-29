# Development Roadmap

## Future

- Nuclear density functional theory interface
- Feynman diagram computation helpers
- Lattice QCD observable interfaces
- Integration with kiran for real-time physics rendering

## Engineering Backlog

### v1.1 — Expanded Wavefunctions and Spectroscopy

- [ ] Hydrogen wavefunctions for n > 3 (general Laguerre polynomial)
- [ ] Angular wavefunctions (spherical harmonics Y_lm)
- [ ] Transition matrix elements (radial integrals)
- [ ] Oscillator strengths from first principles
- [ ] Multi-electron atom approximations (Hartree screening)
- [ ] All spectral series for He+ and Li2+ (Z > 1 validation)

### v1.2 — Nuclear Data Tables

- [ ] AME2020 atomic mass evaluation integration
- [ ] ENDF cross-section data for more isotopes
- [ ] Nuclear charge radii from elastic electron scattering
- [ ] Magnetic dipole and electric quadrupole moments
- [ ] Beta-decay ft values

### v1.3 — Advanced Scattering

- [ ] Partial-wave analysis (phase shifts)
- [ ] Born approximation for screened Coulomb
- [ ] Electron-atom elastic scattering (Mott with form factors)
- [ ] Compton scattering (Klein-Nishina)
- [ ] Pair production cross-sections

### v1.4 — Relativistic Quantum

- [ ] Dirac equation solutions (hydrogen-like)
- [ ] Relativistic corrections beyond first order
- [ ] Hyperfine structure (magnetic dipole interaction)
- [ ] Anomalous magnetic moment corrections
- [ ] Breit interaction for two-electron atoms

### v1.5 — Simulation Support

- [ ] Monte Carlo decay simulation (stochastic chains)
- [ ] Particle transport (energy loss, range, straggling)
- [ ] Detector response functions (Gaussian smearing)
- [ ] Phase space generators for multi-body decays

## Cross-Crate Bridges

- [ ] `bridge.rs` module — primitive-value conversions for cross-crate atomic/nuclear physics
- [ ] **bijli bridge**: electron orbital energy (eV) → photon emission wavelength (nm); nuclear charge → Coulomb field strength
- [ ] **kimiya bridge**: atomic number, electron configuration → valence electrons, electronegativity; isotope mass → molecular weight
- [ ] **prakash bridge**: energy level transitions → spectral line wavelengths and intensities; nuclear spin → hyperfine splitting

## Soorat Integration

- [ ] `integration/soorat.rs` module — feature-gated `soorat-compat`
- [ ] **Atomic orbital visualization**: orbital type (s/p/d/f), quantum numbers, probability density grid for volumetric rendering
- [ ] **Nuclear structure**: nucleon positions (protons/neutrons) for particle rendering
- [ ] **Spectral line data**: wavelength, intensity, element ID for spectral plot rendering
- [ ] **Decay chain**: parent/daughter nuclide graph with half-lives for node-link rendering
