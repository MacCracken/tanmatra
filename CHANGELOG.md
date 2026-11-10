# Changelog

All notable changes to tanmatra will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2.0.2] — 2026-11-09

The roadmap's one unblocked 2.0.x item: the Z > 1 validation. No function changes, so the bundles differ from 2.0.1
only in their version line, and no golden value moves.

### Added
- **The hydrogen-like functions checked against NIST at Z = 2 and 3.** `tests/hydrogenic_nist.tcyr` compares them at
  ⁴He⁺ and ⁷Li²⁺ with the NIST ASD v5.12 He II and Li III energy levels (Yerokhin & Shabaev 2015 for n ≤ 2; Erickson
  1977, scaled by NIST, for n ≥ 3). Each bound is twice the leading term the function leaves out. The new
  `docs/development/reference-validation.md` gives the rule and every derivation, written before the first comparison.
  - `spectral_line_vacuum_nm`, against all 10,026 level pairs with n₁ < n₂ ≤ 20 per ion, within twice the two levels'
    first-order fine structure. The five named series divided by μ match it bit for bit. The worst residual is 0.5004
    of the bound (He⁺) and 0.5005 (Li²⁺).
  - `spectral_line_fine_nm`'s j = l ± ½ splittings for n = 2–11, 55 per ion, within twice the electron's
    anomalous-moment term (α/π)ΔF. Worst: 0.498 and 0.547.
  - `reduced_mass_factor(2, 4)`, against CODATA 2022 m_α/m_e, within 5.8 × 10⁻¹² (twice the dropped 79 eV electron
    binding). Worst: 0.501.

  No defect was found. The 61,948 assertions pass on x86_64 ELF, x86_64 PE and aarch64.
- **The reference data is generated, not typed.** The NIST snapshot is kept byte for byte in `tests/reference/nist_asd/`,
  with its query, citation and SHA-256s. `scripts/reference/nist_levels.py` checks every row (quantum numbers, term,
  parity, completeness for n ≤ 11, the hashes) and writes `tests/nist_levels_data.cyr`. Its `--check` runs in CI.
  Each f64 literal carries its bit pattern, and the test proves the compiler read every one exactly.

### Changed
- **The roadmap lists only open work**, pinned to the release that ships it: 2.0.x patches (two upstream follow-ups
  remain), 2.1.0 atomic structure, 2.2.0 simulation support, the bridges as consumers land, and a 3.0.0 backlog. The
  port plan now records that kimiya's Cyrius port is under way.

## [2.0.1] — 2026-11-09

Two pieces of work that just missed 2.0.0, and finish it:
- **the chrono bridge**, which the plan had already unblocked for the first bridge release (port plan §6);
- **prakash 2.15.2**, which fixes the three bundle warnings that 2.0.0 accepted for the optics profile.

No API is removed or changed, and no golden value moves. The corpus grows by the bridge's 6,588 rows to 191,583, all
within their class on x86_64 ELF, PE and aarch64.

### Added
- **The chrono bridge, the first bridge** (port plan §6), in the new `src/bridge.cyr`: `utc_date_to_tai_offset`,
  `tai_to_utc_seconds` and `utc_to_tai_seconds`, with Rust's names and signatures, over the ported
  `leap_seconds_at`. It is in both bundles.
  - As in Rust, it takes plain values. A `lib/chrono.cyr` datetime supplies them: `dt_year`, `dt_month`, and
    `dt_epoch_secs` + 378,691,200 s for this crate's 1958 epoch.
  - **Verified.** All 6,588 `bridge.tsv` golden rows are 0 ulp on x86_64 ELF, PE and aarch64. `validate.py`
    recomputes the conversions in Python and checks the offsets against `timekeeping.tsv` and the IERS shape.
  - `tests/bridge.tcyr` ports the 4 Rust chrono tests. Three Cyrius-only tests drive the bridge from `lib/chrono.cyr`
    for every month 1972–2035 and check ΔAT against tanmatra's MJD calendar (3,854 assertions). `chrono` joins
    `[test] stdlib`; the library never includes it.
  - Bench `bridge/utc_to_tai_seconds_1000`: 51.8 µs, 12× Rust's `Instant` time (the leap-table scan).
- **Rust-test ledger:** 424 ported, 1 partial, 25 n/a, 48 deferred (was 420 / 52). **Surface parity** handles a
  partly ported file. The bridge has 4 of 4 ported items present, and its other 27 are listed by the consumer each
  waits on. A stale or already-ported entry fails the check.

### Known limitation, filed upstream
- **`lib/chrono.cyr` 6.7.18 gets every date before 1970 wrong, with no error**, so the chrono recipe holds from 1970
  on (`src/bridge.cyr` says so). The bridge's ΔAT is 0 before 1972 anyway. Filed as
  `cyrius/docs/development/issues/2026-11-09-chrono-dates-before-1970-silently-wrong.md`, with a self-proving repro.

### Changed
- **The optics profile pins prakash 2.15.2** (was 2.15.1). The 21 prakash functions that `src/optics.cyr` reaches
  (`spd_new`, `spd_at`, `wavelength_to_rgb` and their callees) are byte-identical across the two releases. All
  1,039 optics golden rows pass on x86_64, PE and aarch64, with the same ulp table, and the optics benches are
  unchanged.
- **hisab is no longer vendored.** prakash 2.15.2 made hisab a test-scope dependency, so `lib/hisab.cyr` (3.6 MB)
  is gone. `cyrius.lock` was re-pinned, and its one commit pin is now prakash; the `sakshi` stdlib leaf stays, with
  the same hash.
- **`dist/tanmatra-optics.deps`**, which `cyrius distlib` derives by compile-verifying the fold, no longer lists
  `sakshi` and `ganita`, which hisab had required. sakshi, the AGNOS logging library, still reaches an optics
  consumer: prakash's own sidecar lists it, its bundle includes it, and the README has consumers list it in
  `[deps] stdlib`, as tanmatra does.

### Fixed
- **`cyrius test` is warning-free in every file again.** prakash 2.15.2 fixed the three bundle warnings that 2.0.0
  accepted in `tests/optics.tcyr` and `tests/optics_parity.tcyr`. tanmatra filed the issue, and prakash archived
  it as fixed.

## [2.0.0] — 2026-11-08

**tanmatra is a Cyrius library.** 2.0.0 ports the Rust crate 1.3.0, frozen in `rust-old/` as the parity oracle, to
Cyrius 6.7.18. Every module is ported except the bridges and the soorat integration, which wait on their consumers
(port plan §6).

### Breaking
- **The language.** The Rust crate is not the 2.0.0 API; it stays at 1.3.0 on crates.io. The Cyrius API follows Rust's names and behavior, with the mapping in the port plan §5:
  - free functions keep their names, an enum's or packed value's methods are `snake(type)_method`, structs are `Tm<Name>` with an `impl`, and globals carry `TM_`;
  - fallible functions return a value-form `Result` with `TM_ERR_*` codes;
  - `Vec` returns fill caller arrays.
- **Not carried over:**
  - serde, except for the five exchange types, whose JSON M7 deferred until a Cyrius consumer needs it;
  - the APIs deprecated in 1.3.0 (`corrected_ft_value`, `FrequencyStandard::quality_factor`);
  - the `std`, `logging` and `full` features. `optics` remains, as `[features] optics` and the `[lib.optics]` profile.

### Verified
- **184,995 golden rows** from the oracle (`tests/golden/`) pass on x86_64 ELF, x86_64 PE and aarch64, in their tolerance classes.
  - Almost all are bit-identical: every exact and arithmetic row, every reaction and scattering row (through bit-exact ports of Rust libm's sin, cos, log, hypot, round and ldexp), and every Einstein A outside one transition.
  - One row is accepted under ADR 0007 (ganita `pow(-1, DBL_MAX)`).
  - The measured class-C bounds and the per-target histograms are in `docs/development/port-parity.md`.
- **The Rust tests.** All 498 Rust tests (476 unit, 20 integration, 2 doc) are accounted for in `docs/development/rust-test-ledger.md`, gated by `scripts/rust_test_ledger.py --check`: 420 ported, 1 partly, 25 not applicable (serde, deprecated), 52 deferred (bridge, soorat, exchange-type JSON).
- **Surface parity:** 254/254 Rust public items have their Cyrius name, and every `rust-old/src` file is ported, deferred or not applicable (`scripts/surface_parity.py`).

The milestone entries below record the port in the order it was done.

### Port — M7 Release 2.0.0 (2026-11-08)
- **The Rust-test ledger** (`scripts/rust_test_ledger.py`, which writes `docs/development/rust-test-ledger.md`; `--check` fails on drift or an unaccounted test).
- **`tests/integration.tcyr`** ports the 16 non-serde tests of `rust-old/tests/integration.rs` and the two doc tests (51 assertions).
- **Surface parity** is the full diff: every file under `rust-old/src` is ported, deferred (bridge, soorat) or not applicable (`lib.rs`'s declarations and prelude), and an unknown file fails it.
- **JSON for the exchange types** stays deferred, as M7 allows: no Cyrius project depends on tanmatra yet, and cyrius issue `2026-11-07-derive-serialize-key-scan-reads-past-nul` is still open.
- **VERSION 2.0.0**; the README describes the Cyrius library; CI runs the ledger check.

### Changed
- **tests**: 8 assertions updated for lints new in clippy 1.99 (`assert_is_empty`, `manual_midpoint`), so the Rust oracle is clean on rustc 1.99.1 before the Cyrius port freezes it.
- **deny.toml**: the license allowlist is trimmed to the four licenses in the dependency graph (MIT, Apache-2.0, GPL-3.0-only, Unicode-3.0).
- **docs**: `docs/development/tanmatra-port-plan.md` revised for Cyrius 6.7.18 (ganita 1.2.15). The 6.6.5 start condition is met.

### Port — M0 steps 1–4 (2026-11-08)
- **`cyrius port`** scaffolded the Cyrius project on 6.7.18. The Rust 1.3.0 crate moved to `rust-old/` as the frozen parity oracle: sources, tests, benches, examples, manifests and `target/`, plus the Rust CI (`rust-old/.github/`), cargo scripts (`rust-old/scripts/`), testing guide and architecture overview (`rust-old/docs/`). The oracle is green from its new location: 476 + 20 + 2 tests, clippy and deny.
- **cyrius.cyml**: `version = "${file:VERSION}"`; stdlib `math`, `ganita` (1.2.15) and `overflow`; `bench` in `[test] stdlib`. `lib/` and `cyrius.lock` are vendored.
- **CLAUDE.md**, **roadmap** and **state.md**: restored or rewritten after the scaffold replaced them, with the Cyrius gates and rules merged in.
- **CI** keeps the scaffold's Cyrius workflow, plus `test --bounds`, `deps --locked` and `audit` steps.
- **ADRs**: `docs/decisions/` → `docs/adr/` (four-digit names). New ADR 0006, the Cyrius port with the Rust oracle, supersedes 0003 (`no_std`).
- **README** and **CONTRIBUTING** describe the port, the Cyrius build and the oracle.
- The scaffold's stub functions (`src/main.cyr`, `tests/tanmatra.{tcyr,bcyr,fcyr}`) have doc comments now, so `cyrius audit` is clean.
- Filed upstream: `cyrius/docs/development/issues/2026-11-08-cyrius-port-clobbers-project-docs-and-omits-bench.md`.

### Port — M0 step 5: golden corpus generator (2026-11-08)
- **`scripts/oracle/`**: `tanmatra-oracle`, a cargo crate that drives the frozen `rust-old/` crate by path. Its lockfile is seeded from the oracle's, and `libm` is pinned `=0.2.16`. It writes `tests/golden/<module>.tsv`, format 1: typed tokens with f64s as IEEE-754 bit patterns, a tolerance class per row and a readable note. The spec is in `scripts/oracle/README.md`.
- Output is deterministic. `--check` fails on any differing, missing or stale file.
- **Golden corpus (M1 plus math reference):**
  - `constants`: all 34;
  - `particle`: all 18 functions over every variant;
  - `relativity`: all 17 functions plus `Add`, over the §3 grid with NaN, ±∞, massless, spacelike and negative-energy edges;
  - `libm`: every libm function tanmatra calls, the reference for the step 8 ulp baselines.
  - 18,453 rows in total.
- **`scripts/oracle/validate.py`** independently re-checks the corpus without the Rust code. It verifies the constants against the source literals, recomputes the 2,962 class-A/B relativity rows bit for bit in Python, and checks sqrt and cbrt.

### Port — M0 step 6: generated constants and tables (2026-11-08)
- **`scripts/oracle/`** also generates the port's data modules from the oracle's API; the binary is now `oracle`, and `--check` covers them. The rules and names are in `scripts/oracle/README.md`.
  - f64s are positional decimal literals, checked to round-trip; there is no hex in a `const`, and ∞ is `(1.0 / 0.0)`.
  - Tables are baked column arrays.
  - Globals carry the `TM_` prefix.
- **`src/constants.cyr`**: all 34 CODATA 2022 constants as `const TM_<NAME>`, carrying the Rust doc citations.
- **`src/particle_data.cyr`**: `TmQuark`, `TmLepton`, `TmBoson` and `TmFundamentalForce`, plus every tabulated PDG 2024 property as a column indexed by enum value; mediators are flattened into start, count and member tables.
- **`scripts/oracle/probe_tables.py`** compiles all 129 generated values through the Cyrius compiler and checks them bit for bit against the golden rows: 0 failures, 0 warnings.

### Port — M0 step 7: parity harness (2026-11-08)
- **`tests/parity.tcyr`** replays every module in `tests/golden/modules.txt`, which the generator now writes. Each row is ok, FAILED, measured (class C without a bound) or pending (not ported).
  - Per function it prints the row counts and an ulp histogram. FAILED rows are shown with their note and both bit patterns, and class-C outliers are reported.
  - The exit code is the number of failing modules.
  - Constants (all 34) and libm sqrt, fabs, floor, trunc and cbrt are bit-identical. sin, cos and tan also match bit for bit; exp, log, sinh and pow are ≤ 1 ulp, and hypot ≤ 2 ulp. All are on x86_64.
- **`tests/ulp.cyr`**: `ulp_bits_match`, `ulp_distance`, `ulp_rel_match` and `ulp_rel_err`, plus the assert wrappers `assert_bits_eq`, `assert_ulp` and `assert_rel`.
- Negative tests (a tampered value, a malformed token, a missing golden file) each fail with exit 1. `cyrius test`, `test --bounds` and `audit` pass with no warnings.
- Filed upstream: `ganita/docs/development/issues/2026-11-08-pow-minus-one-huge-even-exponent.md`. `ganita_f64_pow(-1, y)` returns +∞ or +0 instead of 1 for |y| > 2^64, because `_gn_pow_pos` lacks fdlibm's |x| == 1 guard. tanmatra is unaffected.

### Port — M0 step 8: baselines and class-C bounds; M0 complete (2026-11-08)
- **`docs/development/port-parity.md`**: `libm.tsv` baselines on x86_64 ELF, x86_64 PE (wine) and aarch64 (qemu); macOS is unmeasured.
  - sqrt, fabs, floor, trunc, cbrt, sin, cos and tan are bit-identical to Rust libm on every target.
  - exp is bit-identical on aarch64 and ≤ 1 ulp on x86_64 (x87). log, sinh and pow are ≤ 1 ulp, and hypot ≤ 2.
- **The class-C bound for the libm primitives is 2 ulp**, from both sides' documented < 1 ulp. `tests/parity.tcyr` gates them in ulps (`_par_c_bound`). A check with exp's K = 0 failed exactly the 60 x86_64 rows.
- **ADR 0007** accepts the one `pow(-1, DBL_MAX)` row (ganita's filed defect). The harness lists accepted rows by exact input (`_par_accepted`) and reports them on every run.
- **Hypot.** `invariant_mass` is ill-conditioned for p ≫ m: a 1-ulp change in |p| moves it by more than 16 ulp on 103 of 515 rows. So the port will reproduce Rust libm's hypot as `_tm_hypot`, and the four-momentum golden rows are reclassified from C to B, a tightening. A prototype matched all 600 hypot rows bit for bit on x86_64 and aarch64.
- **Rule for tanmatra's own class-C bounds:** K_F = 2 · S · calls + 2, from F's measured sensitivity, set before F's first comparison. If K_F > 64, make the primitive exact.
- **M0 is complete.** Next is M1.

### Port — M1 Foundation (2026-11-08)
- **`src/numerics.cyr`**: private helpers that reproduce the oracle exactly:
  - `_tm_hypot`, Rust libm's hypot operation for operation, bit-identical on all 600 `libm.tsv` rows on x86_64 ELF, x86_64 PE and aarch64;
  - `_tm_bits` (bit patterns as integers), `_tm_is_inf` and `_tm_nan`, and `_TM_PI` and `_TM_INF`.
- **`src/error.cyr`**: `enum TmError` (`TM_ERR_*`, one per `TanmatraError` variant), `tm_error_name` and `tm_error_message`. The payload and `err_out` design arrives with M2's first fallible API.
- **`src/particle.cyr`**: every `Quark`, `Lepton`, `Boson` and `FundamentalForce` method as `snake(enum)_method(v)`, read from the generated `particle_data.cyr`. An out-of-range value gives NaN, −1 or 0, never an abort. `mediator` is `fundamental_force_mediator_count` plus `fundamental_force_mediator(f, k)`.
- **`src/relativity.cyr`**: `struct TmFourMomentum` with an `impl` block (`TmFourMomentum_new`, `_at_rest` and `_from_mass_and_momentum`; `momentum_magnitude`, `invariant_mass`, `beta`, `gamma` and `kinetic_energy`; `+`), and the nine free functions. Rust's operation order is kept throughout.
- **Parity.** The particle and relativity dispatchers are in `tests/parity.tcyr`. All 6,870 particle and relativity rows are **0 ulp** on all three targets: 17,487 ok, 0 failed, 965 pending (round and ldexp), 1 accepted.
  - `libm::hypot` rows move to `_tm_hypot` and to class A, a tightening.
- **Tests.** 45 Rust tests are ported, in `tests/particle.tcyr` and `tests/relativity.tcyr` (114 assertions), plus a Cyrius-only out-of-range test.
- **Benches.** `tests/tanmatra.bcyr` and `docs/benchmarks-rust-v-cyrius.md`:
  - `lorentz_gamma_1000`: 11.5 µs vs Rust's 4.8 µs;
  - `invariant_mass_1000`: 144 µs vs 21 µs;
  - `particle/properties_all`: 289 ns vs 15 ns.
- **Bundle.** `cyrius.cyml` `[lib]` lists the six modules, and `dist/tanmatra.cyr` is generated by `cyrius distlib`. CI adds `distlib --check` and `scripts/surface_parity.py`, which finds 76/76 ported items.
- **Probe for M2.** A ≤ 8-byte struct in a value-form `Result` can't be bound back safely. So `Nucleus` will be a packed `u64`; filed as `cyrius/docs/development/issues/2026-11-08-small-struct-result-payload-read-as-address.md`.

### Port — M2a Nuclear structure (2026-11-08)
- **`src/nucleus_data.cyr`** (generated): the AME2020 mass excesses (32), charge radii (12), moments (16), superallowed ft and Ft (9), the Mayer–Jensen levels and closures, the empirical proton and neutron level orders, and the SEMF / Myers–Swiatecki / AME2020 u→keV coefficients.
- **Private values.** The coefficients, level orders and `MS_MAGIC` have no public accessor, so `scripts/oracle/src/rust_src.rs` extracts them strictly. The generator then proves each one by recomputing `binding_energy`, `binding_energy_shell_corrected` and `ground_state_spin_parity` over all 14,160 nuclei, and every AME2020 atomic mass, to identical bits. A 1-ulp change to one coefficient makes it refuse to write.
- **`src/nucleus.cyr`**: `Nucleus` as a packed word `z | (a << 32)`:
  - `nucleus_new` returns a value-form Result with `TM_ERR_*`;
  - accessors and presets;
  - SEMF with the Myers–Swiatecki shell term, masses and the radius;
  - the AME2020 and radius lookups as Options, and `nucleus_nuclear_moments(n, &out)`;
  - `TmSuperallowedDecay` with `corrected_ft_seconds`;
  - `superallowed_average_ft` as a pair.
- **`src/shell_model.cyr`**: `ShellLevel` as a packed word with `degeneracy`, `j` and `label`, plus `shell_model_level(k)`, `shell_occupation` into caller buffers, `ground_state_spin_parity` (empirical orders, Brennan–Bernstein rules) and the closures.
- **Golden.** `tests/golden/nucleus.tsv` holds 74,297 rows, with every (Z, A) for 1 ≤ Z ≤ 118, Z ≤ A ≤ 3Z for SEMF, the shell correction and J^π. **All are 0 ulp on x86_64 ELF, x86_64 PE and aarch64.** `Nucleus::atomic_mass_amu` (pow) is class C with K = 2 by the port-parity rule (S = 0). The Myers–Swiatecki term turned out class B.
- **Tests.** `tests/nucleus.tcyr` ports 58 Rust tests (1,765 assertions) plus the surviving assert of `corrected_ft_value_returns_average`.
- **Benches:**
  - `binding_energy_1000`: 135 µs vs 17.7 µs;
  - `binding_energy_shell_corrected_1000`: 1.08 ms vs 133 µs;
  - `ground_state_spin_parity_100`: 13.6 µs vs 0.96 µs, down from 29.2 µs after the level walk was made to read only the 2j column;
  - `shell_occupation_126`: 31 µs vs 13.5 µs.
- **Tooling.**
  - `[lib]` adds `nucleus_data`, `nucleus` and `shell_model`; stdlib adds `tagged` and `result`.
  - Generated `.cyr` files have no trailing space on blank comment lines.
  - The harness takes up to 256 tokens per row.
  - Surface parity is 113/113.

### Port — M2b Decay; M2 complete (2026-11-08)
- **`src/decay_data.cyr`** (generated from `known_isotopes()`): `enum TmDecayMode`, and the 114 NUBASE2020 isotopes as baked columns (Z, A, half-life, mode, isomer flag, excitation energy) plus a NUL-terminated name table with offsets.
- **`src/decay.cyr`:**
  - `decay_constant`, `activity_bq` (Result, `TM_ERR_INVALID_HALF_LIFE`) and `remaining_fraction`;
  - `alpha_decay`, `beta_minus_decay` and `beta_plus_decay` (Result, `TM_ERR_DECAY_NOT_POSSIBLE`);
  - `known_isotope_count` + `known_isotope(k)` (`TmIsotope` by value);
  - `decay_chain` into caller buffers, with the isomeric-transition ground-state rule;
  - `bateman_chain`: the scaling-and-squaring matrix exponential, its scratch `fl_alloc`'d and freed.
- **`src/numerics.cyr`** adds `_tm_ldexp`, Rust libm's scalbn transcribed (bit-identical on all 306 `libm::ldexp` rows, now class A rather than pending), plus `_tm_is_finite` and `_TM_LN2`.
- **Golden.** `tests/golden/decay.tsv` holds 8,511 rows: every isotope, every half-life function with invalid inputs, the α/β grid, chains from every isotope, and the plan's Bateman cases (equal λ, 9-member, the U-238 series at 1 y / 10³ y / 4.463e9 y, the U-235 and Th-232 series). All pass on x86_64 ELF, x86_64 PE and aarch64:
  - class A/B: 0 ulp;
  - `remaining_fraction` (C, K = 4 by the port-parity rule): ≤ 1 ulp;
  - Bateman (D): ≤ 9 ulp on x86_64 (x87 exp) and **0 on aarch64** (the fdlibm polyfill, Rust libm's algorithm).
- **Tests.** `tests/decay.tcyr` ports 24 Rust tests (66 assertions); the 3 serde tests are not carried over.
- **Benches:**
  - `decay_chain_10`: 5.17 µs vs 3.09 µs;
  - `bateman_chain_3`: 5.42 µs vs 1.01 µs;
  - `bateman_u238_series_15`: 270 µs vs 38.2 µs, from 574 µs before the inner product walked pointers;
  - `known_isotopes_all`: 1.93 µs vs 2.36 µs, faster than Rust's allocating `Vec`.
- **Decided.** No error detail travels with a `TM_ERR_*` code: Rust's error strings only restate the caller's input.
- **Tooling.** stdlib adds `freelist`; `[lib]` adds `decay_data` and `decay`; surface parity is 122/122.
- **M2 is complete.** Next is M5 Timekeeping.

### Port — M5 Timekeeping (2026-11-08)
- **`src/timekeeping_data.cyr`** (generated):
  - `enum TmFrequencyStandard` and `enum TmTimeScale`;
  - the CIPM 2025 frequencies and the stabilities;
  - the time-scale constants (`TM_LG_RATE`, `TM_LB_RATE`, `TM_T0_JD`, ...);
  - the 28-entry leap-second table, **derived from the public `leap_seconds_at`** and proven over 1900–2200 × months 0–13;
  - the 13 USNO 1961–1971 segments and `MJD_1972`, extracted and proven over 24,000 quarter-day MJDs.
- **`src/timekeeping.cyr`:**
  - the frequency-standard functions, the TCG/TT and TCB/TDB conversions, and TAI/TT/GPS;
  - `leap_seconds_at` and `tai_minus_utc_seconds_mjd`, with Fliegel & Van Flandern's integer calendar;
  - the relativistic corrections;
  - `struct TmAtomicInstant`: integer seconds and nanoseconds, `+|` saturation and a Euclidean carry. `nanoseconds_since` returns Rust's i128 as (hi, lo), using 128-bit subtract, × 10⁹ through `mulh64`, and carry.
- **`src/numerics.cyr`** adds `_tm_round_half_away` (Rust libm's round; Cyrius `f64_round` rounds halves to even) and `_tm_f64_to_i64_sat` (Rust's saturating `as`).
- **Golden.** `tests/golden/timekeeping.tsv` has 20,703 rows, including exact .5 ns ties found by search. **All are 0 ulp on x86_64 ELF, x86_64 PE and aarch64.** `sagnac_correction_ns` (sin) is class C with K = 6 by the port-parity rule (S = 2). `libm::round` is 659/659, so the golden corpus has no pending rows left: 121,963 ok and 1 ADR-accepted.
- **Tests.** `tests/timekeeping.tcyr` ports 69 Rust tests (193 assertions).
- **Bench.** `atomic_instant_add_seconds_1000`: 33.1 µs vs Rust's 12.9 µs.
- **Generator.** snake() keeps acronyms whole (`TimeScale::TAI` → `TM_TIME_SCALE_TAI`); every earlier generated file is byte-identical. The parity report's function column is wider.

### Port — M3 Atomic (2026-11-08)
- **`src/atomic_data.cyr`** (generated):
  - `enum TmOrbitalType` (with its l, capacity and symbol columns), `enum TmTransitionType`, and `enum TmElectronAffinity`, the bound / unbound / unknown kind;
  - the 118 first ionization energies and 118 electron affinities, through the public API, each labelled with its citation from the Rust source;
  - the Aufbau filling order and the Bethe logarithms (Drake & Swainson 1990), extracted;
  - the 20 NIST ground-state exceptions, **derived from the public `electron_configuration`** as each element's two subshell changes from the Madelung filling.
  - The configurations are proven for Z = 0–120 and u32::MAX, and the Bethe logarithms by recomputing `lamb_shift_nlj_ev` over Z = 0–118 and every (n ≤ 6, l, 2j). A 1-ulp nudge or a dropped exception is refused.
- **`src/atomic_config.cyr`:**
  - the `OrbitalType` properties;
  - `struct TmQuantumNumbers`, validated by `TmQuantumNumbers_new(n, l, ml, ms, out)`;
  - `struct TmOrbitalFilling`, Rust's `OrbitalFilling` with its public fields. A 21-bit packed word failed the golden rows for n = electrons = 2³⁰, so configurations are struct arrays;
  - `electron_configuration` into a caller array, and `format_configuration(_short)`;
  - `ionization_energy_ev`, `electron_affinity` (the kind) and `electron_affinity_ev`.
- **`src/atomic_hydrogenic.cyr`:** Rydberg and reduced-mass lines, first-order fine structure, the five named series (caller buffers), Landé g, Zeeman and Stark, the radial wavefunctions to n = 60, the E1 selection rules, and Einstein A/B from the exact radial dipole integral.
- **`src/atomic_qed.cyr`:** the Lamb shift with Bethe logarithms, vacuum polarization, the Dirac energies, hyperfine splitting, the free and bound g-factors, the strong-field Zeeman shift and Breit.
- **Integer arguments** outside Rust's u32 / i32 give the function's Err, NaN, or `FORBIDDEN`. Where Rust's u32 arithmetic would overflow (`2n − 1` for n ≥ 2³¹), the port is exact in i64.
- **Golden.** `tests/golden/atomic.tsv` has 43,309 rows. **Every class A/B row, and every Einstein A/B row (class D), is 0 ulp on x86_64 ELF, x86_64 PE and aarch64.**
  - Class-C bounds by the port-parity rule, measured by perturbation before the first comparison: `lamb_shift_nlj_ev` and `lamb_shift_ev` K = 10 (S = 4), `radial_wavefunction` K = 6 (S = 2), `radial_probability_density` K = 18 (S = 8).
  - Measured: Lamb ≤ 2 ulp on x86_64 and ≤ 1 on aarch64; R_nl ≤ 2 and R²r² ≤ 7 on x86_64, 0 on aarch64.
  - Dirac turned out class B (sqrt only), a tightening of the plan's C.
  - The whole corpus: 165,272 ok, 0 failed, 0 pending, 1 accepted, on all three targets.
- **`scripts/oracle/validate.py`** independently recomputes the 2,679 class-B line and level rows bit for bit, and checks that every configuration holds Z electrons and is Madelung except at the 20 NIST exceptions.
- **Tests.** `tests/atomic.tcyr` ports 110 of the 114 Rust tests (597 assertions), among them the NIST A-values, the measured Lamb shifts and the radial normalization; the 4 serde round-trips are not carried over.
- **Benches.** 7, mirroring Rust's: 1.9–3.3× Rust for the computed bodies (configurations, radial functions, Einstein A), 6–43× for the call-bound ones (`spectral_line_1000`, `ionization_energy_118`, `lamb_shift_nlj_all_n4`). Two measured optimizations: `#inline` branch-free range checks (`spectral_line_1000` 14.8 → 11.9 µs), and the dipole integral's pow factors and factorial products hoisted out of the double sum, with the same operations and bits (36.1 → 25.4 µs).
- **Surface parity** 41/41 for atomic, 206/206 in total.
- **Generator.** `rust_src.rs` reads `(integer, Enum::Variant)` tables, a pure lookup function's match arms, and per-entry citation comments; a const now ends at its first `;` outside a comment.
- **`src/error.cyr`**'s header no longer promises an `err_out` detail pointer (decided against in M2).

### Port — M6 Optics profile (2026-11-08)
- **`cyrius.cyml`:** an optional `[deps.prakash]` (tag 2.15.1, commit-pinned in cyrius.lock with its hisab and sakshi), `[features] optics = ["prakash"]` (`default = []`, as in Rust), and `[lib.optics]`, the core plus `src/optics_data.cyr` and `src/optics.cyr`. `cyrius distlib --all` writes `dist/tanmatra-optics.cyr`. The core bundle and its sidecar are unchanged and never pull in prakash. A name scan finds no collision between tanmatra's 529 bundle names and prakash's 1,061. `sakshi` joins `[deps]` stdlib, which prakash's bundle needs in scope.
- **`src/optics.cyr`** (the one src file with an include, `lib/prakash.cyr`): `struct TmSpectralLine`; `hydrogen_balmer_series` and `spectral_series_lines` into caller arrays; `lines_to_spd`, a prakash Spd via `spd_new` with every invalid-range guard; `line_to_rgb`, prakash's `wavelength_to_rgb` with its `err_out`.
- **`src/optics_data.cyr`** (generated): the Storey & Hummer Case B Balmer intensities, from the oracle's output, with their upper levels proven against its wavelengths.
- **Golden.** `optics.tsv` has 1,039 rows from the Rust oracle with the Rust prakash 1.2.0, replayed against the Cyrius prakash 2.15.1 on all three targets.
  - `line_to_rgb` (838 rows, every branch edge) and the Balmer series are 0 ulp.
  - `lines_to_spd` is class C with K = 38 (S = 3 by perturbation, 6 exp per sample): ≤ 2 ulp on x86_64, 0 on aarch64.
  - `spectral_series_lines` (class D) has 3 of 2,250 values at 55–95 ulp, about 2 × 10⁻¹⁴ relative. All three are the 49p → 1s Einstein A, where the dipole integral's 49-term sum cancels and ganita's pow is the one inexact primitive. This qualifies M3's "every Einstein row is exact", which held for M3's grid.
- **Harness.** The replay machinery moves to `tests/parity_core.cyr`. `tests/parity.tcyr` keeps the core dispatchers and skips the optics module, so the core suite needs no prakash; `tests/optics_parity.tcyr` replays optics. Rows take up to 1,024 tokens (an SPD row has 804).
- **Tests.** `tests/optics.tcyr` ports the 4 Rust tests (8 assertions); the serde round-trip is not carried over.
- **Benches.** `tests/optics.bcyr`: `balmer_spd_401` 68.1 µs vs Rust's 12.4 µs, and `spectral_series_lyman_10` 14.1 µs vs 7.32 µs (Rust timed with Instant).
- **CI** resolves with `--features optics` and checks the lockfile both ways. Surface parity is 4/4 for optics.
- ⚠ **Known warnings, accepted until prakash fixes them.** prakash 2.15.1's bundle adds three compile warnings to the two optics test files, and tanmatra's code adds none. Two are raw-bit compares that prakash documents as intended (`lib/prakash.cyr:9213–9214`); the third is an unreachable `num_fft` from hisab, which the bundle doesn't include. Filed upstream as `prakash/docs/development/issues/2026-11-08-bundle-consumer-compile-warnings.md`, with two checked warning-free forms of the compare and the options for `num_fft`.

### Port — M4 Reactions & scattering (2026-11-08)
- **Exact sin, cos and log** (`src/numerics.cyr`): `_tm_sin`, `_tm_cos` and `_tm_log` transcribe Rust libm 0.2.16's sin, cos, k_sin, k_cos, rem_pio2 and log operation for operation; arguments past 2²⁰·π/2 use the stdlib's Payne–Hanek reduction. The port-parity rule required this. By perturbation, 12 of the 16 transcendental M4 functions came out past K = 64: the Mott factor at β → 1, the form factor near its zeros, Klein–Nishina's closed form above γ = 1e-2, ξ = 1 + (≈ −1) at large A. `f64_ln` is x87 on x86_64 and a different fdlibm form on aarch64, and `f64_sin` skips Rust's |x| ≤ 9π/4 reductions. A 680,010-input sweep matches Rust libm on every input on x86_64 ELF, x86_64 PE and aarch64; the cost is about 30% over the builtins. `libm.tsv`'s sin, cos and log rows move to the new primitives at class A, and are exact on every target.
- **`src/reaction_data.cyr`** (generated):
  - the seven preset reactions, the eight ENDF/B-VIII.0 thermal cross sections, the U-235 and Pu-239 chain yields, the two nucleosynthesis pathways, `enum TmNucleosynthesisProcess` and `TM_CNO_NEUTRINO_LOSS_MEV`, through the public API; names and step reactions as UTF-8 in one byte table;
  - the private resonance-integral table and scattering's ζ(3), extracted and proven.
- **`src/reaction.cyr`:**
  - `struct TmNuclearReaction`, with products and leptons as word arrays (`r.product(k)`, `r.lepton(k)`), the projectile 0 for Rust's None, and `nuclear_charge_change` / `baryon_number_change`;
  - the seven presets;
  - `q_value` over two buffers, its sums folded from −0.0 as Rust's `Sum for f64` does;
  - the Coulomb barrier, geometric cross section, Breit–Wigner and the resonance integrals;
  - the thermal cross sections and yields as a count and an accessor returning a struct by value;
  - the pathways, with `p.step(k)`;
  - neutron moderation.
- **`src/scattering.cyr`:** Rutherford and Mott, partial waves (phase shifts as a buffer and a count), Born for a screened Coulomb potential, the Thomas–Fermi length, Mott with the uniform form factor (series below |qR| = 0.5), Compton and Klein–Nishina (series below γ = 1e-2), and Maximon's pair production (both branches).
- **Golden.** `reaction.tsv` (7,965 rows) and `scattering.tsv` (10,704) cover every preset, every branch point on both sides (with the adjacent doubles), and θ from 0 to 2π plus ±∞ and NaN. **All 18,669 rows are 0 ulp on x86_64 ELF, x86_64 PE and aarch64.** Every function except `rutherford_total_above_angle` is class B; that one uses ganita's tan, class C with K = 10 (S = 4), and is also 0 ulp. The whole corpus: 183,956 ok, 0 failed, 0 pending, 1 accepted.
- **Tests.** `tests/reaction.tcyr` ports 39 of the 43 Rust reaction tests (133 assertions; the 4 serde round-trips are not carried over), and `tests/scattering.tcyr` all 53 scattering tests (770 assertions), among them the Hubbell–Gimm–Øverbø pair-production value and the Simpson-rule checks of the Rutherford and Klein–Nishina totals.
- **Benches.** `rutherford_1000` 62.2 µs vs Rust's 6.55 µs; `pair_production_1000` 101.5 µs vs 16.6 µs. A sin call is most of the gap.
- **Surface parity** 26/26 for reaction and 18/18 for scattering.
- **Generator.** `cyr.rs` writes i64 tables (nucleus words, enum names).
- **Bundle.** `[lib]` adds reaction_data, reaction and scattering. `dist/tanmatra.deps` now lists `math`: `_tm_rem_pio2` calls the stdlib's Payne–Hanek reduction for huge arguments.

### Port — `ShellLevel` as a struct (2026-11-08)
- **`src/shell_model.cyr`**: `ShellLevel` is `struct TmShellLevel { n_shell; l; two_j; }` with `impl` `new`, `degeneracy`, `j` and `label`, in place of M2's 21-bit packed word and its `shell_level_*` functions. Rust's fields are public u32s, so a hand-built level past 2²¹ is valid there, and the word silently gave it the wrong label, j and degeneracy (11 of the new golden rows failed against it). `shell_model_level(k)` returns the struct by value, and `shell_occupation` fills a caller `TmShellLevel[32]`.
- **Golden.** `nucleus.tsv` gains 5 hand-built probes (15 rows) with fields at 2²¹, 2³⁰ and u32::MAX: 74,312 rows. Everything passes on x86_64 ELF, x86_64 PE and aarch64: 165,287 ok, 0 failed, 1 accepted.
- **Benches.** `ground_state_spin_parity_100` 13.6 → 11.5 µs (the walk uses table indices, no level values); `shell_occupation_126` 31.3 → 11.6 µs, now 0.86× Rust.
- **Tests and harness.** The nucleus tests, parity dispatchers and bench use the struct; their level buffers are typed `TmShellLevel[32]` (768 bytes): the old 256-byte buffers held 32 words but would hold only 10 structs. Surface parity moves `ShellLevel` to the struct rules: 37/37.

## [1.3.0] — 2026-09-16

Repairs every finding of the 2026-09-16 math and data audit
(`docs/audit/2026-09-16-math-audit.md`). Every data table was regenerated from
its primary source and every formula fix is pinned by a reference-value test.
This is the final Rust feature release and the reference implementation for the
Cyrius port, which will ship as 2.0.0 (`docs/development/tanmatra-port-plan.md`).
It contains API breaks, listed below; the version number follows the port plan
rather than a Rust major bump.

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
