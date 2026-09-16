# tanmatra Rust → Cyrius port plan

> **Priority work.** This is the plan of record for porting tanmatra to Cyrius.
> Rust **1.3.0** is the final Rust feature release and the **reference implementation**;
> the Cyrius port ships as **2.0.0**.
>
> ⏸ **Start condition: wait for the Cyrius 6.6.5 release.** tanmatra needs `f64_cbrt`,
> which first shipped in ganita 1.2.6. Cyrius 6.6.4 still bundles ganita 1.2.5, so the
> port begins on 6.6.5, the release that bundles ganita ≥ 1.2.6. The pre-port assessment
> found no other blocker; the cube-root gap it exposed is resolved upstream.

| Thing | Value | Source |
|---|---|---|
| Rust source (oracle) | tanmatra **1.3.0**, 11,418 lines in `src/` | `Cargo.toml`, `VERSION` |
| Rust health at port start | 476 unit + 20 integration + 2 doc tests; fmt, clippy `-D warnings`, rustdoc `-D warnings`, `cargo audit` and `cargo deny` all clean; 19 criterion benches | `docs/guides/testing.md` |
| Correctness baseline | Every finding of the 2026-09-16 audit repaired and pinned by reference-value tests | [`../audit/2026-09-16-math-audit.md`](../audit/2026-09-16-math-audit.md) |
| Cyrius toolchain | **6.6.5** — the port starts on this release (pin in `cyrius.cyml`). 6.6.4 bundles ganita 1.2.5, which lacks `f64_cbrt` | cyrius releases |
| ganita | **1.2.6** (first release with `f64_cbrt`), reached through the 6.6.5 stdlib bundle | `/home/macro/Repos/ganita` |
| prakash | **2.4.0**, Cyrius (ported) | `/home/macro/Repos/prakash` |
| bijli, kimiya, soorat | **Rust**, not yet ported (kimiya 1.1.1, soorat 1.0.0) | respective repos |
| Target | tanmatra **2.0.0**, Cyrius | this document |

Porting model to copy: **prakash**. It is a science library with the same shape as
tanmatra: flat modules, heavy f64 math, generated constant tables, and a
Rust-vs-Cyrius benchmark record.
- `/home/macro/Repos/prakash/cyrius.cyml` — `[lib]` bundle order
- `src/spectral_cie.cyr` — generated hex tables
- `docs/benchmarks-rust-v-cyrius.md`

Language and toolchain references:
- `/home/macro/Repos/cyrius/docs/guides/cyrius-guide.md`
- `docs/stdlib-reference.md`
- `docs/api-surface.snapshot`
- `/home/macro/Repos/szal/docs/development/port-plan.md` §1 — Rust→Cyrius cheat-sheet

---

## 1. Scope

### In scope for 2.0.0

| Rust module | Code lines | Public fns | Tests | Cyrius module(s) |
|---|---|---|---|---|
| `error.rs` | 36 | — | 2 | `src/error.cyr` |
| `constants.rs` | 134 | (constants) | — | `src/constants.cyr` |
| `particle.rs` | 294 | 18 | 28 | `src/particle.cyr` |
| `relativity.rs` | 232 | 17 | 22 | `src/relativity.cyr` |
| `nucleus.rs` | 940 | 37 | 63 | `src/nucleus_data.cyr`, `src/nucleus.cyr`, `src/shell_model.cyr` |
| `decay.rs` | 570 | 9 | 27 | `src/decay_data.cyr`, `src/decay.cyr` |
| `atomic.rs` | 1,684 | 41 | 114 | `src/atomic_data.cyr`, `src/atomic_config.cyr`, `src/atomic_hydrogenic.cyr`, `src/atomic_qed.cyr` |
| `reaction.rs` | 780 | 25 | 43 | `src/reaction_data.cyr`, `src/reaction.cyr` |
| `scattering.rs` | 543 | 18 | 53 | `src/scattering.cyr` |
| `timekeeping.rs` | 747 | 35 | 74 | `src/timekeeping.cyr` |
| `optics.rs` (feature `optics`) | 172 | 4 | 5 | `src/optics.cyr`, opt-in `[lib.optics]` profile against **prakash 2.4.0** |

prakash is already in Cyrius, so the `optics` integration is ported with parity:
- `lines_to_spd` → prakash `spd_new`
- `line_to_rgb` → prakash `wavelength_to_rgb(wavelength_nm, err_out)`
- `hydrogen_balmer_series` and `spectral_series_lines` are pure tanmatra math

As with prakash's `ai` module, `optics.cyr` ships as an opt-in bundle profile so a
consumer that doesn't need prakash doesn't pull it in.

### Deferred until the consumers exist in Cyrius

| Rust surface | Why deferred | Unblocked when |
|---|---|---|
| `bridge.rs` — bijli bridges (`energy_to_wavelength_nm`, `nuclear_charge_field`) | bijli is not in Cyrius | bijli ports |
| `bridge.rs` — kimiya bridges (`atomic_number_to_valence`, `neutron_count`, `mass_to_binding_deficit_mev`, `half_life_to_decay_constant`) | kimiya is Rust (1.1.1) | kimiya ports |
| `bridge.rs` — prakash-labelled helpers (`transitions_to_wavelengths`, `nuclear_spin_to_hyperfine_scale`) | kept with the bridge module for one coherent later port | bridge phase |
| `bridge.rs` — jyotish, chrono, hisab-mimamsa, falak bridges and `SimulationClock` (kiran/joshua) | consumers not in Cyrius | each consumer ports |
| `bridge.rs` — bhava `TimeContext` | bhava is Rust (2.0.0) | bhava ports |
| `integration/soorat.rs` (feature `soorat-compat`) | soorat is Rust (1.0.0) | soorat ports |

Every deferred function is a thin wrapper over in-scope physics, for example
`tai_to_tt` over timekeeping and `mass_to_binding_deficit_mev` over constants.
Porting the physics first means each bridge later is a small, parity-testable
file. `rust-old/` keeps the Rust bridge code and its tests as the reference until
then.

### Not carried over

- **Deprecated in 1.3.0:** `nucleus::corrected_ft_value` (ignored its argument) and
  `FrequencyStandard::quality_factor` (unsourced values). Their replacements
  `SuperallowedDecay::corrected_ft_seconds`, `superallowed_average_ft` and
  `quality_factor_for_linewidth` are in scope.
- **Rust-only mechanics:** `no_std`/`alloc` gating, `Send + Sync` assertions, and
  cargo features, which become `[lib]` profiles.
- **serde:** the Rust suite's serde round-trip tests do not port as-is. JSON
  exchange is ported only for types a Cyrius consumer actually exchanges, via
  `#derive(Serialize)`/bayan — see M7.

---

## 2. M0 — Initialize with `cyrius port`

`cyrius port` is the initializer. Do not hand-create the project structure; if the
tool misses something, fix the tool.

```sh
cd /home/macro/Repos/tanmatra
cyrius port --dry-run .     # review: moves src/, Cargo.toml, tests/, benches/ → rust-old/
cyrius port .               # scaffolds src/main.cyr, cyrius.cyml (pins the installed toolchain — must be 6.6.5), CI, port doc templates
```

After the scaffold:

1. **`cyrius.cyml`**
   - `[package].cyrius = "6.6.5"`.
   - `[deps] stdlib`: `string, fmt, alloc, io, vec, str, syscalls, assert, bench, math`, plus `ganita`.
   - `[lib] modules` in dependency order: error → constants → particle → relativity → nucleus_data → nucleus → shell_model → decay_data → decay → atomic_* → reaction_* → scattering → timekeeping. `optics` goes in `[lib.optics]`.
2. **ganita cube root — resolved: wait for 6.6.5.** tanmatra calls `cbrt` in the binding energy, the nuclear radius, Thomas–Fermi screening and form factors. `f64_cbrt` first appears in ganita **1.2.6** (filed from the tanmatra assessment as `ganita/docs/development/issues/archived/2026-09-16-f64-cbrt-missing.md`). Cyrius 6.6.4 bundles ganita 1.2.5, so M0 starts only on **Cyrius 6.6.5**, which bundles ganita ≥ 1.2.6. There is no separate `cyrius deps` pin for ganita. At M0, confirm `lib/ganita.cyr` in the 6.6.5 snapshot reports version ≥ 1.2.6 and that `f64_cbrt` is present.
3. **Golden-value generator in `rust-old/`** — a Rust example, e.g. `rust-old/examples/golden.rs`, that writes `tests/golden/<module>.tsv`: one line per call, with inputs and outputs as **IEEE-754 hex bit patterns** plus the tolerance class (§4). It covers every in-scope public function over the input grids listed per milestone. It is the parity oracle and is regenerated only from Rust 1.3.0.
4. **Table generator in `rust-old/`** — writes the Cyrius data modules (`*_data.cyr`) with every f64 as exact hex, following prakash `spectral_cie.cyr`. Tables are generated, never re-typed; the audit found 17/32 AME rows and ~30/118 EA rows hand-transcribed wrong.
5. **Parity harness** — `tests/parity.tcyr` loads `tests/golden/*.tsv` and asserts each row within its class. Counts are printed per module and the exit code is the number of failing modules (ecosystem repro convention).
6. **Baseline measurements** — before porting any formula, run the golden inputs through Cyrius `f64_exp`, `f64_ln`, `f64_sin`, `f64_cos`, `f64_pow`, `f64_cbrt` and `f64_sinh` against Rust `libm` and record the ulp distributions in `docs/development/port-parity.md`. These fix the tolerance of class C (§4) from evidence, not guesswork.

**Gate:** `cyrius build --strict` of the empty bundle, `cyrius test` running the parity harness (0 modules yet), both generators committed.

---

## 3. Milestones

Each milestone ends with the same gates:
- `cyrius build --strict`
- `cyrius test`, with the module's golden rows plus its ported reference-value tests
- `cyrius bench`, with the module's benches recorded against the Rust numbers
- `cyrius vet`, `cyrius lint`, `cyrius fmt --check`
- CHANGELOG and `state.md` updated

| Milestone | Modules | Key work | Golden grid (minimum) |
|---|---|---|---|
| **M1 Foundation** | error, constants, particle, relativity | Constants as hex (CODATA 2022). `TanmatraError` → error enum. `FourMomentum` as a 4-f64 heap struct with accessors. Enum methods (`Quark`, `Lepton`, `Boson`, `FundamentalForce`) as `fn quark_mass_mev(q)` switch tables | every enum variant; β ∈ {0, ±0.5, ±(1−1e-12), ±1, ±1.5}; four-momenta from 1 keV to 10 GeV |
| **M2 Nuclear** | nucleus_data, nucleus, shell_model, decay_data, decay | Generated AME2020, radii, moments, superallowed and 114-isotope tables. SEMF + Myers–Swiatecki (`f64_cbrt`). Mayer–Jensen and empirical p/n orders; Brennan–Bernstein coupling. `decay_chain` with the IT ground-state rule. **Bateman scaling-and-squaring matrix exponential** | all (Z, A) with 1 ≤ Z ≤ 118, Z ≤ A ≤ 3Z (binding, shell-corrected, J^π); all isotopes; the U-238/U-235/Th-232 chains; Bateman: equal λ, 9-member long-lived chain, U-238 series at 1 y / 10³ y / 4.463e9 y |
| **M3 Atomic** | atomic_data, atomic_config, atomic_hydrogenic, atomic_qed | Generated IE (118) and EA (118, 3-state) tables. Configurations with the 20 NIST exceptions. Rydberg and reduced-mass lines. Dirac and first-order fine structure. Radial functions (Laguerre coefficients and factorial-ratio normalisation) for n ≤ 60. Exact radial dipole integrals → Einstein A/B. Lamb shift with Bethe logarithms, VP, hyperfine, g-factors, Breit | Z = 1..118 (configs, IE, EA); n1 < n2 ≤ 12, Z ∈ {1, 2, 3, 26, 92}; R_nl for all n ≤ 10, l < n, 50 radii; A-values for every allowed n_u ≤ 10 → n_l; all (n, l, j) with n ≤ 4 |
| **M4 Reactions & scattering** | reaction_data, reaction, scattering | ENDF/B-VIII.0 thermal cross sections, resonance integrals and yields (generated). `NuclearReaction` layout with optional projectile, neutron counts and lepton vec. Presets. Rutherford, Mott (p·β), Born, form factor (series branch), Klein–Nishina (series branch), Maximon pair production (both branches) | all presets; energies spanning each branch point (KN γ = 1e-2; form factor qR = 0.5; pair k = 4) on both sides; θ ∈ {1e-3, 0.1, π/2, π} |
| **M5 Timekeeping** | timekeeping | CIPM 2025 standards. Leap-second table plus the 1961–1971 USNO segments and MJD→calendar. `AtomicInstant` in integer nanoseconds. Geoid-referenced satellite correction. Sagnac (loop and one-way). TCG/TT and TCB/TDB conversions | leap lookups 1961–2030; `AtomicInstant` add/diff at epochs 0, 2.1e9 s and ±9.2e18 s edges; J2000 and T0 conversions |
| **M6 Optics profile** | optics (prakash 2.4.0) | `lines_to_spd` → `spd_new`; the Storey & Hummer Balmer series; A-value series intensities; `line_to_rgb` → `wavelength_to_rgb` | Balmer SPD at 0.1 and 1 nm FWHM; invalid-range guards |
| **M7 Release 2.0.0** | all | JSON for exchange types (`Nucleus`, `Isotope`, `FourMomentum`, `AtomicInstant`, `NuclearReaction`) if a consumer needs them. `cyrius distlib` bundle. `docs/benchmarks-rust-v-cyrius.md`. `cyrius audit`. VERSION 2.0.0 | full golden corpus green |

Suggested order: M1 → M2 → M5 → M3 → M4 → M6 → M7. Timekeeping is self-contained
and exercises the integer-time paths early. Atomic is the largest module.

---

## 4. Parity policy

Every golden row carries one tolerance class:

| Class | Tolerance | Applies to |
|---|---|---|
| **A — exact** | bit-identical, or equal integers/strings | constants and table lookups; configurations; J^π; decay chains; leap seconds; `AtomicInstant`; magic numbers; enum properties |
| **B — arithmetic** | ≤ 2 ulp | formulas using only + − × ÷ √ and correctly rounded `f64_cbrt`: SEMF, Rutherford, Klein–Nishina (series branch), relativity, Q-value arithmetic |
| **C — transcendental** | relative ≤ bound measured in M0 step 6 (expected ~1e-14) | anything through exp/ln/pow/trig: shell term, Dirac, Lamb shift, Mott, pair production, form factor, radial functions |
| **D — iterative** | relative ≤ 1e-12 | Bateman matrix exponential, radial dipole integrals, numerically integrated references |

Rules:
- **Never widen a class to make a row pass.** A class-C failure beyond the M0 bound is a defect: in the port, in the Cyrius math function (file it upstream as with `f64_cbrt`), or in the Rust oracle (fix Rust, bump 1.3.x, regenerate).
- **The physics reference tests travel with the code.** The audit's reference-value tests (NIST A-values, measured Lamb shifts, AME2020 binding energies, IERS offsets, ENDF values) are ported as Cyrius assertions in addition to golden parity: parity proves "same as Rust", the reference tests prove "right".
- **Branch points are tested on both sides:** the Klein–Nishina series switch, form-factor series switch, pair-production k = 4 join, Bateman degenerate λ, and the pre/post-1972 TAI−UTC switch.

---

## 5. Rust → Cyrius mapping (tanmatra specifics)

### Math functions

| Rust (`libm` / core) | Cyrius | Note |
|---|---|---|
| `sqrt`, `exp`, `log`, `sin`, `cos`, `atan`, `floor`, `ceil`, `fabs` | `f64_sqrt`, `f64_exp`, `f64_ln`, `f64_sin`, `f64_cos`, `f64_atan`, `f64_floor`, `f64_ceil`, `f64_abs` | builtins |
| `cbrt` | `f64_cbrt` | ganita ≥ 1.2.6 via the Cyrius 6.6.5 stdlib (§2 step 2); correctly rounded, matches musl |
| `pow` | `f64_pow` | exact path for integer exponents; `x^(5/3)` stays `x·cbrt(x)²` as in Rust |
| `sinh`, `hypot`, `trunc` | `f64_sinh`, `f64_hypot`, `f64_trunc` | ganita / stdlib |
| `round` | `f64_round` | confirm ties-away-from-zero matches `libm::round` in the M0 baseline |
| `tan` | local `sin/cos` (hisab `f64_util.cyr` precedent) | used only in Rutherford total |
| `ldexp` | multiply by an exact power of two built from bits | Bateman scaling |
| `f64::powi` | repeated multiplication | Lamb shift `(Zα)⁴` |
| `f64::midpoint` | `(a + b) / 2` | j + ½ |
| `i128` (`AtomicInstant::nanoseconds_since`) | split seconds/nanos arithmetic with checked ops (`+?`, `-?`) | no i128 in Cyrius |

### Data model

| Rust | Cyrius |
|---|---|
| `Nucleus { z: u32, a: u32 }` (Copy, validated) | packed value `(z << 16) \| a` in one i64 with `nuc_new/z/a/n` accessors: no allocation in hot loops; validation in `nuc_new` returns an error code |
| `#[non_exhaustive] enum` (DecayMode, OrbitalType, TimeScale, FrequencyStandard, …) | `enum` with explicit numeric values; property functions via `switch` |
| `ElectronAffinity::{Bound(f64), Unbound, Unknown}` | tag + value pair (`return (tag, ev)`) or value-form `Result`; not a heap tagged union in hot lookups |
| `Result<T, TanmatraError>` | value-form `Result` (zero-alloc since cyrius 6.6.0), with a `TanmatraError` code enum plus a detail string only on error. Match hisab 3.x, which moved to `Result<T,E>` at 3.0.0 |
| Structs with identity (`Isotope`, `NuclearReaction`, `ShellLevel`, `SuperallowedDecay`, `ThermalCrossSection`, `FissionYield`, `FourMomentum`, `AtomicInstant`) | heap blob + `#derive(accessors)` + layout comment |
| `Vec<T>` returns (configurations, shell occupation, series, yields, chains, Bateman populations) | `lib/vec.cyr` of pointers, or f64 buffers with explicit length; hot APIs take a caller buffer |
| `known_isotopes()` (allocates 114 structs per call) | lazily built static table + indexed accessors, allocated once, as prakash CIE tables |
| `String` names and labels | `Str` via `str_builder` |

### Cyrius constraints that change tanmatra APIs

- **More than 6 arguments.** `breit_wigner_cross_section` (7) and `born_screened_coulomb_with_masses` (7) exceed the 6 register arguments; take a params struct pointer.
- **Negative literals.** None exist; write `(0 - x)`. Negative table entries (mass excesses, magnetic moments, quadrupoles) are generated as hex bit patterns anyway.
- **`var buf[N]` is static.** Never return one. Bateman and radial-integral scratch use `alloc` or a caller buffer; a long-running consumer loop must not leak through the bump allocator, so scratch goes in an arena.
- **Globals.** At most 4,096 globals with non-literal initializers per compilation unit (raised at v6.3.41; integer-literal initializers and enum members don't count). Keep f64 constants in `fn` accessors returning hex, and tables in lazily initialized buffers, as prakash does. Handle `alloc_reset()` invalidating lazily built tables; prakash `spectral_cie.cyr` documents the hazard and its fix.
- **Unsigned semantics.** u32 has no Cyrius equivalent, so validate non-negative Z, A, n, l and j explicitly at each API boundary.

---

## 6. After 2.0.0

1. **Bridges** (`src/bridge.cyr`) — ported per consumer as each lands in Cyrius:
   - kimiya, bijli, bhava;
   - jyotish, chrono, hisab-mimamsa, falak;
   - kiran/joshua (`SimulationClock`).

   Parity comes from the Rust bridge tests in `rust-old/`.
2. **soorat integration** (`src/soorat.cyr`) when soorat ports.
3. **Removing `rust-old/`** happens only after every deferred bridge and integration is ported or explicitly dropped, and the golden corpus and generators are archived (prakash removed `rust-old/` in 2.2.3 and keeps it recoverable from git).
4. **Roadmap items continue in Cyrius:**
   - v1.1 remainder: Y_lm, oscillator strengths, Hartree screening, He⁺/Li²⁺ validation;
   - v1.6: simulation support.

---

## 7. Tracking

| Milestone | Status | Golden rows | Reference tests | Benches | Notes |
|---|---|---|---|---|---|
| M0 Initialize | ⏸ | — | — | — | waiting for Cyrius 6.6.5 (bundles ganita ≥ 1.2.6 with `f64_cbrt`) |
| M1 Foundation | ☐ | | | | |
| M2 Nuclear | ☐ | | | | |
| M5 Timekeeping | ☐ | | | | |
| M3 Atomic | ☐ | | | | |
| M4 Reactions & scattering | ☐ | | | | |
| M6 Optics profile | ☐ | | | | |
| M7 Release 2.0.0 | ☐ | | | | |
| Bridges (deferred) | ⏸ | | | | waits on bijli, kimiya, bhava, jyotish, falak, hisab-mimamsa, kiran/joshua |
| soorat integration (deferred) | ⏸ | | | | waits on soorat |
