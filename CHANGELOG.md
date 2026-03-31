# Changelog

All notable changes to tanmatra will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0]

### Added

#### Nuclear Data Tables (v1.2)
- **nucleus**: AME2020 atomic mass evaluation — 32 nuclides with experimental mass excess values (keV), `experimental_mass_excess_kev()`, `experimental_atomic_mass_amu()`
- **nucleus**: Nuclear charge radii — 12 nuclides from Angeli & Marinova 2013, `charge_radius_fm()`
- **nucleus**: Nuclear electromagnetic moments — `NuclearMoments` struct, 16 nuclides with magnetic dipole (μ_N) and electric quadrupole (barn) from Stone 2005/2019
- **nucleus**: Superallowed beta-decay ft values — `SuperallowedDecay` struct, 9 transitions from Hardy & Towner 2020, `corrected_ft_value()` (Ft = 3072.27 s)
- **reaction**: Resonance integrals — 12 isotopes from ENDF/B-VIII.0, `resonance_integral_barns()`

#### Advanced Scattering (v1.3)
- **scattering**: Partial-wave analysis — `legendre_polynomial()`, `partial_wave_cross_section()`, `partial_wave_differential()`
- **scattering**: Born approximation for screened Coulomb — `born_screened_coulomb()`, `thomas_fermi_screening_fm()`
- **scattering**: Electron-atom elastic scattering — `mott_electron_differential()`, `nuclear_form_factor_uniform()`, `mott_electron_with_form_factor()`
- **scattering**: Compton scattering — Klein-Nishina formula: `compton_energy_ratio()`, `klein_nishina_differential()`, `klein_nishina_total()`
- **scattering**: Pair production — Bethe-Heitler: `pair_production_cross_section()`

#### Relativistic Quantum (v1.4)
- **atomic**: Dirac equation solutions — `dirac_energy_mev()`, `dirac_binding_energy_ev()` for hydrogen-like atoms
- **atomic**: Relativistic corrections — `relativistic_correction_ev()` captures all orders beyond non-relativistic
- **atomic**: Hyperfine structure — `hyperfine_splitting_ev()` for s-states with nuclear g-factor
- **atomic**: Anomalous magnetic moment — `electron_g_factor()`, `bound_electron_g_factor()`, `anomalous_zeeman_splitting_ev()`
- **atomic**: Breit interaction — `breit_interaction_ev()` for helium-like ion ground states

#### Frequency Standards & Atomic Time (v1.5)
- **timekeeping** (new module): `FrequencyStandard` enum — Cs-133, Rb-87, H-maser, Sr-87 optical, Yb-171 optical with frequencies, wavelengths, quality factors, Allan deviation
- **timekeeping**: `TimeScale` enum — TAI, UTC, TT, GPS, TCB, TCG with conversion functions
- **timekeeping**: Leap second table — 28 entries (1972–2017) from IERS Bulletin C, `leap_seconds_at()`
- **timekeeping**: `AtomicInstant` — TAI-referenced instant (i64 s + u32 ns, epoch 1958-01-01)
- **timekeeping**: Relativistic clock corrections — `gravitational_redshift()`, `schwarzschild_clock_correction_us_per_day()`, `second_order_doppler_shift()`, `sagnac_correction_ns()`
- **bridge**: Jyotish bridge `tai_to_tt()` for exact TAI-TT conversion
- **bridge**: Chrono bridge — `utc_date_to_tai_offset()`, `tai_to_utc_seconds()`, `utc_to_tai_seconds()` for leap-second-aware TAI↔UTC
- **bridge**: Hisab-mimamsa bridge — `gravitational_time_dilation()`, `gravitational_time_offset_s()` for Schwarzschild satellite clock corrections
- **bridge**: Falak bridge — `tai_seconds_to_jd_tt()`, `jd_tt_to_tai_seconds()`, `tai_seconds_to_mjd_tt()` for atomic time ↔ Julian Date
- **bridge**: Kiran/Joshua bridge — `SimulationClock` type with time-scale multiplier, pause/resume/fast-forward
- **bridge**: Bhava bridge — `TimeContext` enum (RealTime/Simulated/Paused) for circadian/rhythm time awareness

#### Constants
- `CLASSICAL_ELECTRON_RADIUS_FM` (CODATA 2022)
- `PROTON_G_FACTOR` (CODATA 2022)
- `ELECTRON_ANOMALOUS_MOMENT` (CODATA 2022)
- `STANDARD_GRAVITY`, `GM_EARTH`, `EARTH_ROTATION_RAD_S`

### Changed
- Test count: 243 → 429 unit tests + 20 integration + 2 doctests

## [1.1.0]

### Added
- **bridge** — cross-crate primitive-value bridges for bijli (energy to wavelength, nuclear charge field), kimiya (valence electrons, neutron count, binding energy deficit, decay constant), prakash (transition energies to wavelengths, nuclear spin to hyperfine)
- **integration/soorat** — feature-gated `soorat-compat` module with visualization data structures: `OrbitalVisualization` (hydrogen-like probability density slice), `NuclearStructure` (proton/neutron positions on golden spiral shells), `SpectralLineData` (wavelength/intensity lines), `DecayChainVisualization` (nuclide/transition graph)

### Updated
- prakash 1.1.0 -> 1.2.0, zerocopy 0.8.47 -> 0.8.48

## [1.0.0] - 2026-03-28

### Added

- **constants**: CODATA 2022 physical constants (17 constants including Bohr magneton, reduced Planck constant in MeV*s)
- **particle**: Standard Model with quarks (6 flavors, PDG 2024 masses), leptons (6 types), gauge bosons (photon, gluon, W, Z, Higgs), four fundamental forces, decay widths and lifetimes (W, Z, Higgs, muon, tau), `is_stable()` method
- **nucleus**: Bethe-Weizsacker binding energy with Strutinsky shell correction, nuclear shell model (Mayer-Jensen, 32 levels through 184), ground-state spin-parity, shell occupation, shell closure functions, magic numbers, nuclear radii, preset nuclei (H-1, He-4, C-12, Fe-56, U-235, U-238)
- **decay**: 114 known isotopes with NNDC/NUBASE 2020 half-lives, complete Th-232, U-235, and U-238 natural decay chains, nuclear isomers (Ta-180m, Hf-178m2, Tc-99m, Am-242m, Pa-234m), isomeric transition decay mode, Bateman equations for sequential decay chain populations, decay constant, remaining fraction, activity
- **atomic**: Electron configurations (Aufbau + 22 NIST exceptions through Z=118), NIST ionization energies (Z=1-118), electron affinities (Z=1-118), hydrogen radial wavefunctions (n=1-3), radial probability densities, Einstein A/B coefficients, electric dipole selection rules, fine-structure energy levels, named spectral series (Lyman, Balmer, Paschen, Brackett, Pfund), Rydberg spectral lines, Zeeman effect (anomalous, Lande g-factor), Stark effect (linear, hydrogen), Lamb shift, vacuum polarization
- **reaction**: Q-value computation, Coulomb barrier estimation, 7 preset reactions (DT, DD x2, pp, U-235 fission, CNO, triple-alpha), cross-sections (geometric, Breit-Wigner resonance, thermal neutron database with 8 key isotopes), neutron moderation (lethargy, collisions to thermalize, moderating ratio), fission product yield distributions (U-235, Pu-239 from ENDF/B-VIII.0), nucleosynthesis pathways (s-process, r-process)
- **relativity**: Four-momentum with `Add` trait, Lorentz factor, velocity addition, relativistic energy-momentum relations, de Broglie wavelength, invariant mass calculations, velocity conversions
- **scattering**: Rutherford differential and total cross-sections, Mott scattering correction, distance of closest approach, Sommerfeld parameter
- **optics**: Feature-gated prakash integration (`optics` feature) with `SpectralLine` type, spectral line to SPD conversion with Gaussian profiles, Balmer series generator, spectral series generator, wavelength-to-RGB color mapping
- **error**: `TanmatraError` with 6 variants, full serde + thiserror support
- Full serde support on all 19 public types with roundtrip tests
- `no_std` + `alloc` support
- `#![forbid(unsafe_code)]`
- Send + Sync assertions on all key types
- 11 criterion benchmarks
- 242 tests with 91.7% coverage
- 4 runnable examples (basic, nuclear, spectral, relativity)
- Architecture decision records (5 ADRs)
- Threat model documentation
- Testing guide with benchmark results
- CI/CD pipelines (check, test, coverage, MSRV, release)
