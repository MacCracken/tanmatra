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

### v1.2 — Nuclear Data Tables (shipped in 1.2.0)

- [x] AME2020 atomic mass evaluation integration
- [x] ENDF cross-section data for more isotopes
- [x] Nuclear charge radii from elastic electron scattering
- [x] Magnetic dipole and electric quadrupole moments
- [x] Beta-decay ft values

### v1.3 — Advanced Scattering (shipped in 1.2.0)

- [x] Partial-wave analysis (phase shifts)
- [x] Born approximation for screened Coulomb
- [x] Electron-atom elastic scattering (Mott with form factors)
- [x] Compton scattering (Klein-Nishina)
- [x] Pair production cross-sections

### v1.4 — Relativistic Quantum (shipped in 1.2.0)

- [x] Dirac equation solutions (hydrogen-like)
- [x] Relativistic corrections beyond first order
- [x] Hyperfine structure (magnetic dipole interaction)
- [x] Anomalous magnetic moment corrections
- [x] Breit interaction for two-electron atoms

### v1.5 — Frequency Standards & Atomic Time (shipped in 1.2.0)

Atomic transitions define the SI second. tanmatra already models transitions (Einstein A/B, selection rules, Rydberg) — this version adds the specific transitions and time scale conversions that define physical timekeeping.

#### Frequency Standards

- [x] `FrequencyStandard` enum: Cesium133, Rubidium87, HydrogenMaser, StrontiumOptical, YtterbiumOptical
- [x] Cesium-133 hyperfine transition: 9,192,631,770 Hz (SI second definition, CGPM 1967)
- [x] Strontium-87 optical lattice: 429,228,004,229,873.2 Hz (BIPM 2021 secondary representation)
- [x] Ytterbium-171 optical lattice: 518,295,836,590,863.6 Hz (BIPM 2021)
- [x] Rubidium-87 hyperfine: 6,834,682,610.904 Hz
- [x] Hydrogen maser 1420.405 MHz (21 cm line)
- [x] `transition_frequency()`, `transition_wavelength()`, `quality_factor()` for each standard
- [x] Fractional stability (Allan deviation) characterization per standard type

#### Atomic Time Scales

- [x] `TimeScale` enum: TAI, UTC, TT, GPS, TCB, TCG
- [x] TAI↔UTC conversion with leap second table (IERS Bulletin C, maintained as const data)
- [x] TAI↔TT: fixed offset TT = TAI + 32.184s (IAU 1991)
- [x] TAI↔GPS: fixed offset GPS = TAI − 19s
- [x] TCB↔TCG↔TT relativistic coordinate time conversions (IAU 2000)
- [x] `AtomicInstant` type: TAI-referenced, sub-nanosecond precision (i64 seconds + u32 nanos)
- [ ] Conversions to/from `chrono::DateTime<Utc>` via leap second table

#### Relativistic Clock Corrections

- [x] Gravitational redshift at altitude: Δf/f = −gΔh/c² (first order)
- [x] Full Schwarzschild correction for satellite clocks (GPS: +45.850 μs/day gravitational, −7.214 μs/day velocity)
- [x] Second-order Doppler shift for moving clocks
- [x] Sagnac correction for rotating reference frames (Earth surface)
- [ ] Bridge to hisab-mimamsa: `gravitational_time_dilation()` for exact corrections

#### Cross-Crate Integration (Bridges)

- [x] **jyotish bridge**: `tai_to_tt()` replaces jyotish's polynomial Delta T approximation with exact TAI↔TT for modern epochs (post-1958)
- [x] **chrono bridge**: `utc_date_to_tai_offset()`, `tai_to_utc_seconds()`, `utc_to_tai_seconds()` — leap-second-aware TAI↔UTC conversion for chrono consumers
- [x] **hisab-mimamsa bridge**: `gravitational_time_dilation()`, `gravitational_time_offset_s()` — Schwarzschild clock corrections (GPS +38.6 μs/day net)
- [x] **falak bridge**: `tai_seconds_to_jd_tt()`, `jd_tt_to_tai_seconds()`, `tai_seconds_to_mjd_tt()` — atomic time ↔ Julian Date conversion for ephemeris
- [x] **kiran/joshua bridge**: `SimulationClock` type with time-scale multiplier, pause/resume/fast-forward — game world time anchored to physical TAI
- [x] **bhava bridge**: `TimeContext` enum (RealTime/Simulated/Paused) consuming `SimulationClock` — circadian/rhythm/growth modules distinguish wall-clock from simulation time

### v1.6 — Simulation Support

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
