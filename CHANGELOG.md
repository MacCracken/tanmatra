# Changelog

All notable changes to tanmatra will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] — 2026-09-16

Repairs every finding of the 2026-09-16 math and data audit
(`docs/audit/2026-09-16-math-audit.md`). Every data table was regenerated from
its primary source and every formula fix is pinned by a reference-value test.
Major version: `NuclearReaction` gained fields and `projectile` became
`Option<Nucleus>`; several functions now return physically different values.

### Breaking
- **reaction**: `NuclearReaction.projectile` is `Option<Nucleus>` (`None` for neutron-induced reactions); new fields `neutrons_in`, `neutrons_out`, `leptons_out`; `nuclear_charge_change()`, `baryon_number_change()`. `u235_fission` no longer uses H-1 as the neutron.
- **bridge**: `nuclear_spin_to_hyperfine_scale` returns I + 1/2 (was 2I + 1); `atomic_number_to_valence` is derived from the ground-state configuration (was 2 for every Z ≥ 19).
- **nucleus**: `binding_energy()` uses AME2020-fitted coefficients; `binding_energy_shell_corrected()` uses the Myers–Swiatecki shell term; `atomic_mass_amu()` now includes electrons (use `nuclear_mass_amu()` for the nucleus).
- **atomic**: `electron_affinity_ev` table replaced; `lamb_shift_ev`, `vacuum_polarization_ev`, `einstein_a_coefficient`, `breit_interaction_ev`, `radial_wavefunction` return corrected values (see Fixed).

### Fixed — formulas
- **scattering**: `distance_of_closest_approach` was 2× too small; `rutherford_total_above_angle` was 4× too small; `pair_production_cross_section` was negative from threshold to 3.42 MeV (now Maximon's full-range Bethe–Heitler formulas); `mott_electron_differential` uses pc·β from the kinetic energy; Klein–Nishina total and uniform form factor use series where the closed forms cancel; `mott_correction_factor` documented as the spin-½ Mott factor.
- **atomic**: R₃₁ was 6× too large — radial functions are now exact for any n ≤ 60; Einstein A is the exact hydrogenic dipole rate (was off by 0.025–21×, scaled as Z¹⁰); vacuum polarization coefficient 4/15 (was 1/3); Lamb shift from the one-loop Zα expansion with Bethe logarithms (was a Z⁴ scaling with an invented p-state factor); hyperfine splitting includes (I + ½); Breit interaction is the leading α²Z³/4 E_h Breit–Pauli term (was an unsourced Z⁴ formula); fine structure uses the CODATA Rydberg energy (was 13.6 eV); Ds and Rg configurations follow the relativistic predictions (6d⁸7s², 6d⁹7s²); invalid quantum numbers are rejected.
- **nucleus**: odd-odd nuclei no longer get half-integer spins (Brennan–Bernstein coupling rules); empirical proton/neutron level orders (odd-A agreement 38% → 49%); `corrected_ft_value` deprecated (it ignored its argument).
- **decay**: `bateman_chain` is a non-negative scaling-and-squaring matrix exponential — exact for equal decay constants, independent of the time unit, accurate for trace daughters; `decay_chain` no longer loops on isomeric transitions; Bi-212 follows its β⁻ branch.
- **timekeeping**: `AtomicInstant::add_seconds` is exact to 1 ns (was lossy through f64, −73 ns at today's epoch); serde normalizes nanoseconds; Sagnac documented as the closed-loop formula plus new one-way `sagnac_one_way_ns`; GPS correction references the geoid potential (38.575 µs/day); `fractional_stability` values sourced.
- **relativity**: `lorentz_gamma` symmetric in β; precision-preserving forms for γ ↔ β, invariant mass and kinetic energy; `velocity_addition` returns NaN only when undefined.
- **optics**: `lines_to_spd` validates its range (no panic on step 0); Balmer intensities are Storey & Hummer (1995) Case B; series intensities from exact A-values.
- **soorat**: orbital slices use exact R_nl and |Y_l0|² (radial nodes, angular shape); protons and neutrons interleaved.

### Fixed — data
- **constants**: CODATA 2022 values (were CODATA 2018; ħ in MeV·s was 2014, g_p 2010).
- **particle**: PDG 2024 quark masses and Higgs mass (were PDG 2022/2023).
- **nucleus**: AME2020 mass excesses regenerated (17 of 32 were wrong by up to 7.3 keV); F-19 quadrupole moment 0 (spin ½); Ni-58 radius; moments from INDC(NDS)-0794 (2019) and -0833 (2021); superallowed ft and Ft from Hardy & Towner 2020 (were the 2009/2015 surveys).
- **decay**: half-lives regenerated from NUBASE2020 with 1 y = 365.2422 d (62 of 114 differed); Ta-180m has no observed decay (was a lower limit stored as a half-life); Cu-64 dominant mode EC.
- **reaction**: ENDF/B-VIII.0 thermal cross sections with free-atom scattering (H-1 was the 82 b bound value); capture, fission and (n,α) resonance integrals separated (U-235 and B-10 "capture" values were fission and (n,α)); ENDF/B-VIII.0 chain yields; Q-values from AME2020 (U-235 channel 173.28 MeV, CNO 26.731 MeV); r-process third peak Pt-195.
- **atomic**: ionization energies from NIST ASD (Tc was 0.16 eV off; superheavy values from Smits et al. 2023); electron affinities from AHH99 plus later measurements, with an explicit `ElectronAffinity` type distinguishing bound, unbound and unknown.
- **timekeeping**: CIPM 2025 recommended frequencies; TAI−UTC for 1961–1971; leap-second table validity date.

### Added
- **atomic**: `spectral_line_vacuum_nm`, `reduced_mass_factor`, `lamb_shift_nlj_ev`, `hyperfine_splitting_spin_ev`, `electron_affinity`, `ElectronAffinity`.
- **nucleus**: `nuclear_mass_amu`, `SuperallowedDecay::corrected_ft_seconds`, `superallowed_average_ft`.
- **reaction**: `fission_resonance_integral_barns`, `n_alpha_resonance_integral_barns`, `CNO_NEUTRINO_LOSS_MEV`.
- **scattering**: `born_screened_coulomb_with_masses`.
- **timekeeping**: `tcg_minus_tt_seconds`, `tt_to_tcg_jd`, `tcg_to_tt_jd`, `tcb_to_tdb_jd`, `tdb_to_tcb_jd`, `tcb_minus_tcg_secular_seconds`, `LB_RATE`, `TDB0_S`, `T0_JD`, `tai_minus_utc_seconds_mjd`, `LEAP_SECOND_TABLE_VALID_UNTIL`, `sagnac_one_way_ns`, `time_dilation_shift_exact`, `AtomicInstant::add_nanoseconds`, `AtomicInstant::nanoseconds_since`, `FrequencyStandard::quality_factor_for_linewidth`.
- **constants**: `HBAR_C_MEV_FM`, `HC_EV_NM`, `RYDBERG_EV`, `ELECTRON_MASS_U`, `PROTON_MASS_U`, `NEUTRON_MASS_U`, `HYDROGEN_ATOM_MASS_U`, `NUCLEAR_MAGNETON_EV_T`, `DEUTERON_G_FACTOR`, `ATOMIC_UNIT_TIME_S`, `COULOMB_K_SI`.
- 8 benchmarks (19 total); reference-value tests (476 unit + 20 integration + 2 doc).

### Deprecated
- `nucleus::corrected_ft_value` (ignored its argument).
- `FrequencyStandard::quality_factor` (unsourced order-of-magnitude values).

### Security
- crossbeam-epoch 0.9.18 → 0.9.21 (dev dependency via criterion; RUSTSEC-2026-0204).

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
