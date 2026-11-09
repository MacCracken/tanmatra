# tanmatra Rust → Cyrius port plan

> **Priority work.** This is the plan of record for porting tanmatra to Cyrius.
> Rust **1.3.0** is the final Rust feature release and the **reference implementation**;
> the Cyrius port ships as **2.0.0**.
>
> ▶ **Start on Cyrius 6.7.18.** The original start condition was ganita ≥ 1.2.6, for
> `f64_cbrt`. It is met: 6.7.18 bundles ganita **1.2.15**, and its `f64_cbrt` is correctly
> rounded. This revision (2026-11-08) re-checks the plan against the installed 6.7.18
> toolchain. Every language claim in §5 marked *probed* was compiled and run on it.

| Thing | Value | Source |
|---|---|---|
| Rust source (oracle) | tanmatra **1.3.0**, 11,418 lines in `src/`, unchanged since 2026-09-16 | `Cargo.toml`, `VERSION` |
| Rust health, re-run 2026-11-08 on rustc 1.99.1 | ✅ all green after M0 step 0: 476 unit + 20 integration + 2 doc tests; fmt, clippy `-D warnings`, rustdoc `-D warnings`, `cargo audit` and `cargo deny`, with zero warnings. Step 0 fixed 8 test-only lints new in clippy 1.99 and trimmed 7 unused license allowances | local run |
| Correctness baseline | Every finding of the 2026-09-16 audit is repaired and pinned by reference-value tests | [`../audit/2026-09-16-math-audit.md`](../audit/2026-09-16-math-audit.md) |
| Cyrius toolchain | **6.7.18** (installed; `cyrius port` pins it) | `~/.cyrius/current` |
| ganita | **1.2.15**, bundled in the 6.7.18 stdlib. There is no separate pin | `~/.cyrius/versions/6.7.18/lib/ganita.cyr:2` |
| prakash | **2.15.1**, Cyrius, pinned to 6.7.18. `spd_new` and `wavelength_to_rgb` keep their 2.4.0 signatures | `/home/macro/Repos/prakash` |
| Deferred consumers | Still in Rust and untouched since spring: kimiya 1.1.1, soorat 1.0.0, bhava 2.0.0. bijli, jyotish, falak, hisab-mimamsa, kiran and joshua have no local repo | respective repos |
| Target | tanmatra **2.0.0**, Cyrius | this document |

## 0. What changed since the 6.6.5 draft

The scope, the milestones and the parity policy still stand. Cyrius 6.6.6 – 6.7.18 changed
the mechanics underneath them:

- **The start condition is met.** ganita 1.2.15 has a correctly rounded `f64_cbrt`. Rust's
  `libm` 0.2.16 `cbrt` is the CORE-MATH port, also correctly rounded, so cube roots are
  bit-identical.
- **`cyrius port` overwrites project files.** It replaces `CLAUDE.md`, the roadmap,
  `.gitignore` and both workflows, sets `version = "0.1.0"`, and leaves `bench` out of
  `[test] stdlib`. M0 step 2 restores the files, and the bugs are filed upstream.
- **Several §5 constraints are gone:**
  - functions with 7 or more arguments work (stack args), so no params structs;
  - negative literals exist;
  - there is no globals cap;
  - array locals live on the stack.
- **New traps for f64 code.** These four are all *probed*. `(0 - x)` on an f64 is now wrong.
  A hex literal in a `const` is an integer. `f64_round` rounds ties to even. `+?` aborts the
  process.
- **Better data representations exist:**
  - generated tables become **baked typed column arrays**, not lazy hex accessors;
  - `FourMomentum`, `AtomicInstant` and `Nucleus` become **struct values**;
  - `ElectronAffinity` was to be a **`: stack` payload enum**; it can't nest inside a value-form
    `Result`, so it is the kind plus `electron_affinity_ev` (M3, §5);
  - errors propagate with **`?`**.
- **`optics` becomes a real feature:** `[features]` plus `optional = true` on the prakash
  dependency, as well as the `[lib.optics]` bundle profile.
- **Gates changed:**
  - `--strict` is dropped (it has been a no-op since 6.3.2);
  - `cyrius test` must be **warning-free**, because the f64 type diagnostics appear only
    there and not in `cyrius lint`;
  - new gates: `--bounds`, `distlib --check`, `deps --locked` and a surface-parity check.
- **Class C is measured per target.** On x86_64, `f64_exp` and `f64_ln` use x87, so their bits
  differ between targets.
- **ranga joins prakash as a porting model**, for the differential oracle. The szal
  cheat-sheet is retired: it uses the 6.1.33 nomenclature.

## Porting models and references

- **prakash 2.15.1**, `/home/macro/Repos/prakash`. It has tanmatra's shape: a science library
  with flat modules, heavy f64 math, generated constant tables and a Rust-vs-Cyrius
  benchmark record.
  - Copy the `cyrius.cyml` `[lib]` / `[lib.ai]` layout. A profile repeats the whole module
    list; the order is load-bearing, because no `src/` file has an `include`.
  - Copy `docs/benchmarks-rust-v-cyrius.md`.
  - **Don't** copy its lazily built hex tables in `src/spectral_cie.cyr`, or its `alloc_reset`
    workaround. Both predate baked typed arrays (6.6.16).
  - Its `rust-old/` was removed in 2.2.3.
- **ranga 2.3.2**, `/home/macro/Repos/ranga`, pinned to 6.7.18 and keeping `rust-old/`. It is
  the differential-oracle model.
  - `scripts/oracle/`: a Rust driver builds against `rust-old/` by path, and a Cyrius probe
    runs beside it. Both print f64 bit patterns, and `compare.py` reports ulp distances.
    Accepted divergences are listed in its `ACCEPTED` table, each with the ADR that decided
    it.
  - `docs/development/parity-rust-v-cyrius.md` and `docs/development/port-mechanics.md`.
  - ⚠ port-mechanics §3 predates two fixes. Decimal literals have been correctly rounded
    since 6.6.10, and `f64_cbrt` exists now.
- **mela 1.5.2**: `docs/development/rust-parity.md`, a per-module ledger of Rust tests (each
  covered or not applicable). It is the model for the M7 test ledger.
- **Language and toolchain**, in `/home/macro/Repos/cyrius`:
  - `docs/guides/cyrius-guide.md`, `docs/stdlib-reference.md`, `docs/api-surface.snapshot`;
  - `docs/development/ecosystem-migration.md`: read the 6.7.x entries at every pin bump;
  - `cyrius help manifest` for the manifest keys.

---

## 1. Scope

### In scope for 2.0.0

| Rust module | Code lines | Public fns | Tests | Cyrius module(s) |
|---|---|---|---|---|
| `error.rs` | 36 | — | 2 | `src/error.cyr` |
| `constants.rs` | 134 | (constants) | — | `src/constants.cyr` |
| `particle.rs` | 294 | 18 | 28 | `src/particle_data.cyr` (generated), `src/particle.cyr` |
| `relativity.rs` | 232 | 17 | 22 | `src/relativity.cyr` |
| `nucleus.rs` | 940 | 37 | 63 | `src/nucleus_data.cyr`, `src/nucleus.cyr`, `src/shell_model.cyr` |
| `decay.rs` | 570 | 9 | 27 | `src/decay_data.cyr`, `src/decay.cyr` |
| `atomic.rs` | 1,684 | 41 | 114 | `src/atomic_data.cyr`, `src/atomic_config.cyr`, `src/atomic_hydrogenic.cyr`, `src/atomic_qed.cyr` |
| `reaction.rs` | 780 | 25 | 43 | `src/reaction_data.cyr`, `src/reaction.cyr` |
| `scattering.rs` | 543 | 18 | 53 | `src/scattering.cyr` |
| `timekeeping.rs` | 747 | 35 | 74 | `src/timekeeping.cyr` |
| `optics.rs` (feature `optics`) | 172 | 4 | 5 | `src/optics.cyr`, behind the `optics` feature and the `[lib.optics]` profile, against **prakash 2.15.1** |

The counts were re-verified 2026-11-08. They are unchanged.

prakash is already in Cyrius, so the `optics` integration is ported with parity:
- `lines_to_spd` → prakash `spd_new(start_nm, step_nm, values, len)`, which borrows `values`
  and returns 0 if allocation fails;
- `line_to_rgb` → prakash `wavelength_to_rgb(wavelength_nm, err_out)`, which sets
  `PK_ERR_WAVELENGTH_RANGE` outside 380 – 780 nm;
- `hydrogen_balmer_series` and `spectral_series_lines` are pure tanmatra math.

The default bundle must not pull in prakash, as with prakash's own `ai` module. §2 step 3
shows the feature wiring.

### Deferred until the consumers exist in Cyrius

| Rust surface | Why deferred | Unblocked when |
|---|---|---|
| `bridge.rs`: bijli bridges (`energy_to_wavelength_nm`, `nuclear_charge_field`) | bijli is not in Cyrius (no repo yet) | bijli ports |
| `bridge.rs`: kimiya bridges (`atomic_number_to_valence`, `neutron_count`, `mass_to_binding_deficit_mev`, `half_life_to_decay_constant`) | kimiya is Rust (1.1.1) | kimiya ports |
| `bridge.rs`: prakash-labelled helpers (`transitions_to_wavelengths`, `nuclear_spin_to_hyperfine_scale`) | kept with the bridge module for one coherent later port | bridge phase |
| `bridge.rs`: jyotish, hisab-mimamsa and falak bridges, and `SimulationClock` (kiran/joshua) | consumers not in Cyrius | each consumer ports |
| `bridge.rs`: chrono bridge (`utc_date_to_tai_offset`, `utc_to_tai_seconds`, `tai_to_utc_seconds`) | Rust targets the `chrono` crate | **unblocked**: the Cyrius stdlib has `lib/chrono.cyr`. Port it in the first bridge release |
| `bridge.rs`: bhava `TimeContext` | bhava is Rust (2.0.0) | bhava ports |
| `integration/soorat.rs` (feature `soorat-compat`) | soorat is Rust (1.0.0) | soorat ports |

Every deferred function is a thin wrapper over in-scope physics; for example, `tai_to_tt`
wraps timekeeping and `mass_to_binding_deficit_mev` wraps constants. Porting the physics
first makes each later bridge a small file that can be parity-tested. Until then,
`rust-old/` keeps the Rust bridge code and its tests as the reference.

### Not carried over

- **Deprecated in 1.3.0.** `nucleus::corrected_ft_value` ignored its argument, and
  `FrequencyStandard::quality_factor` used unsourced values. Their replacements are in scope:
  `SuperallowedDecay::corrected_ft_seconds`, `superallowed_average_ft` and
  `quality_factor_for_linewidth`.
- **Rust-only mechanics:** `no_std`/`alloc` gating and `Send + Sync` assertions. Cargo features
  map to `[features]` plus `[lib.<profile>]`.
- **serde.** The Rust suite's serde round-trip tests don't port as they are. JSON exchange is
  ported only for types a Cyrius consumer actually exchanges; see M7.

---

## 2. M0 — Initialize with `cyrius port`

`cyrius port` is the initializer. Don't hand-create the project structure. Where the tool
gets something wrong, repair the output here and file the defect upstream (step 2).

0. ✅ **Freeze a green oracle (Rust, before the move). Done 2026-11-08.**
   - Fixed the 8 test-only lints new in clippy 1.99:
     - 7 × `assert_is_empty`, now `assert_eq!`/`assert_ne!` against `[]`, so a failure prints
       the value;
     - 1 × `manual_midpoint` in `form_factor_small_qr_no_cancellation`. The qR = 0.5 switch
       point is now a named local, which gives the same values bit for bit.
   - Trimmed `deny.toml`'s license allowlist to the four licenses actually in the graph: MIT,
     Apache-2.0, GPL-3.0-only and Unicode-3.0. The other seven were never matched.
   - Only tests and lint config changed, so `VERSION` stays **1.3.0**.
   - The full P(-1) cleanliness check is clean on rustc 1.99.1.
1. ✅ **Scaffold. Done 2026-11-08.**
   ```sh
   cd /home/macro/Repos/tanmatra
   cyrius port --dry-run .
   cyrius port .
   ```
   - **Moved into `rust-old/`:** `Cargo.toml`, `Cargo.lock`, `src/`, `tests/`, `benches/`,
     `examples/`, `Makefile`, `codecov.yml`, `deny.toml` and `target/`. `--dry-run` lists
     none of these by name.
   - **Written:** `src/main.cyr`, `cyrius.cyml` (pinned to `cyrius = "6.7.18"`, with no
     `[lib]`), `.gitattributes` (`lib/** -text`), `tests/tanmatra.{tcyr,bcyr,fcyr}`,
     `docs/{adr,architecture,guides,examples}/` and `docs/development/state.md`. `lib/` and
     `cyrius.lock` come from the first `cyrius deps` (step 3).
   - **Left alone:** `CHANGELOG.md`, `README.md`, `docs/audit/`, this plan and `scripts/`.
2. ✅ **Repair what the scaffold overwrites. Done 2026-11-08.**
   - **`CLAUDE.md` is replaced.** It was merged from `git show HEAD:CLAUDE.md` and the
     scaffold's Cyrius rules, keeping the project identity, consumers, process and DO NOT
     list. The Rust key principles were translated:

     | Rust principle | Cyrius equivalent |
     |---|---|
     | `#[non_exhaustive]` | append-only enum values that are never renumbered |
     | `#[must_use]` | `#must_use`, plus `#pure` where it applies |
     | serde on every type | JSON only for exchange types (M7) |
     | zero unwrap/panic | no `+?`, `-?`, `*?`, `assert_fatal` or `result_unwrap` in `src/` |

   - **`docs/development/roadmap.md` is replaced** (123 → 36 lines). It was restored from
     git. `docs/development/state.md` was rewritten for the real version (1.3.0, not 0.1.0).
   - **`.github/workflows/{ci,release}.yml` are replaced.**
     - The scaffold's Cyrius CI is kept, extended with `test --bounds`, `deps --locked` and
       `audit` steps.
     - The Rust workflows are kept for reference in `rust-old/.github/workflows/`; GitHub
       doesn't run them there.
     - The golden corpus is regenerated locally, so CI has no Rust job.
     - `.gitignore` is appended to, not replaced: `/rust-old/target/` and `/build/`.
   - **Other Rust-only material moved to `rust-old/`:**
     - `scripts/bench-history.sh` and `scripts/version-bump.sh` → `rust-old/scripts/`;
     - `docs/guides/testing.md` and `docs/architecture/overview.md` → `rust-old/docs/`.
   - **ADRs.** `docs/decisions/001`–`005` → `docs/adr/0001`–`0005`, the scaffold's four-digit
     convention. ADR 0006 records the port decision and supersedes 0003 (`no_std`).
   - **Filed upstream:**
     `cyrius/docs/development/issues/2026-11-08-cyrius-port-clobbers-project-docs-and-omits-bench.md`.
     - It overwrites `CLAUDE.md`, the roadmap, `state.md` and the workflows.
     - It writes `version = "0.1.0"` despite `VERSION`, and its own `release.yml` then fails the
       first tag.
     - It leaves `bench` out of `[test] stdlib`.
     - Its stub files have six undocumented functions, so `cyrius audit` fails (exit 1) until
       they get doc comments. tanmatra's stubs are documented now.
     - `--dry-run` doesn't list what it moves or overwrites.
3. ✅ **`cyrius.cyml`. Done 2026-11-08** for everything except `[lib]` (M1) and
   `[deps.prakash]` / `[features]` / `[lib.optics]` (M6). `cyrius deps` vendors `lib/` and
   writes `cyrius.lock`. Both are committed, as in prakash. The shape below comes from
   `cyrius help manifest` and prakash's manifest. Confirm the `[features]` / `optional` keys
   when M6 adds the profile; prakash doesn't use them, so this is new ground.
   ```toml
   [package]
   name = "tanmatra"
   version = "${file:VERSION}"
   license = "GPL-3.0-only"
   language = "cyrius"
   cyrius = "6.7.18"

   [deps]
   stdlib = ["syscalls", "string", "alloc", "str", "fmt", "io", "vec", "math", "ganita", "overflow"]

   [deps.prakash]
   git = "https://github.com/MacCracken/prakash.git"
   tag = "2.15.1"
   modules = ["dist/prakash.cyr"]
   optional = true

   [features]
   optics = ["prakash"]

   [test]
   stdlib = ["assert", "bench"]

   [lib]
   modules = ["src/error.cyr", ...]   # dependency order (below)

   [lib.optics]
   modules = [ <the whole [lib] list>, "src/optics.cyr" ]
   ```
   - `math` is required for `f64_sin` and `f64_cos` on every target, and for
     exp/ln/atan on aarch64.
   - `overflow` provides the saturating `+|`, `-|` and `*|` operators.
   - `[lib]` order: error → constants → particle_data → particle → relativity → nucleus_data → nucleus →
     shell_model → decay_data → decay → atomic_data → atomic_config → atomic_hydrogenic →
     atomic_qed → reaction_data → reaction → scattering → timekeeping.
   - `[lib] modules = []` is invalid, so `[lib]` gets its first entry, `src/error.cyr`, in M1.
4. ✅ **ganita check. Done 2026-11-08.** `lib/ganita.cyr` reports `1.2.15`, and
   `ganita_f64_cbrt`, `_pow`, `_sinh`, `_tan` and `_hypot` are present. Use the canonical `ganita_f64_*` names throughout. The short
   aliases `f64_sinh`, `f64_pow`, `f64_atan2` and `f64_hypot` are marked *Deprecated alias …
   migration window only*.
5. ◐ **Golden-value generator in `scripts/oracle/`. Built 2026-11-08; module coverage
   grows per milestone.** The format spec and the rules are in
   [`scripts/oracle/README.md`](../../scripts/oracle/README.md).
   - **Crate.** `tanmatra-oracle` depends on `rust-old/` by path, with `features =
     ["optics"]`; it never builds inside `rust-old/` and never modifies it. ranga's
     `scripts/oracle/rust/` works the same way.
     - Its `Cargo.lock` was seeded from `rust-old/Cargo.lock`, and `libm` is pinned
       `=0.2.16`, so every shared dependency is the oracle's version.
     - clippy `-D warnings` is clean.
   - **Output.** It writes `tests/golden/<module>.tsv`, format 1:
     - one row per call: `fn<TAB>class<TAB>inputs<TAB>outputs<TAB>note`;
     - typed tokens, with every f64 as its **IEEE-754 bit pattern** (`f:3ff0000000000000`);
     - structs flattened by field, sequences as `v:<len>` followed by their elements;
     - a readable note per row, which the harness ignores.
   - **Deterministic.** Two runs are byte-identical, and `--check` exits 1 on any differing,
     missing or stale file. Grids use only the pinned `libm` crate and seeded xorshift.
   - **Independently checked.** `scripts/oracle/validate.py` (Python, no packages):
     - re-parses every row;
     - checks `constants.tsv` against the Rust source literals;
     - recomputes the 2,962 arithmetic-only class-A/B relativity rows with Python floats in Rust's operation
       order: bit-identical, 0 failures;
     - checks `sqrt` against `math.sqrt` and `cbrt` on perfect cubes.
   - **Coverage so far:**
     - **M1:** `constants` (all 34), `particle` (all 18 functions, every variant: 102 rows),
       `relativity` (all 17 functions plus `Add`, on the §3 grid widened with NaN, ±∞,
       massless, spacelike and negative-energy cases: 6,768 rows);
     - **`libm`** (11,549 rows): every `libm` function tanmatra calls, over special values,
       the ranges tanmatra uses and seeded wide sweeps. It is the reference for step 8.
   - **The remaining modules are added at the start of their milestone** (M2 nucleus, shell
     model and decay; M5 timekeeping; M3 atomic; M4 reaction and scattering; M6 optics), each
     with its §3 grid. The corpus is 2.8 MB so far.
   - It is the parity oracle and is regenerated only from Rust 1.3.0.
6. ◐ **Table and constant generator in `scripts/oracle/`. Built 2026-11-08, in the same
   crate and binary as step 5; data modules grow per milestone.** Rules and names are in
   [`scripts/oracle/README.md`](../../scripts/oracle/README.md). Tables are generated, never
   re-typed: the audit found 17 of 32 AME rows and about 30 of 118 EA rows hand-transcribed
   wrong.
   - **Sources.** Values come from the oracle's public API. Primary data files are not
     re-parsed: the audit already verified the oracle against them, and the golden corpus
     pins the port to the oracle. The Rust source is read only for `///` doc comments, which
     carry the citations into the Cyrius file.
   - **Generated so far:**
     - **`src/constants.cyr`:** all 34 CODATA 2022 constants as `const TM_<NAME>`.
     - **`src/particle_data.cyr`:**
       - the enums `TmQuark`, `TmLepton`, `TmBoson` and `TmFundamentalForce`, with values in
         Rust declaration order;
       - every tabulated property as a column indexed by enum value, e.g.
         `TM_QUARK_MASS_MEV[TM_QUARK_UP]`;
       - `FundamentalForce::mediator` flattened into start, count and member tables.
       - Properties that compute ħ over a tabulated value stay computed, in class B.
   - **Literals.**
     - f64s are decimal literals written positionally, using the shortest round-trip
       digits. Each one is checked in Rust to parse back to the oracle's bits.
     - The lexer has **no exponent notation**: `1.5e-3` is a parse error.
     - Decimal literals have been correctly rounded since 6.6.10. For example,
       `0.0000000000000000001602176634` gives `0x3c07a4da290c1653`, the bits of Python's
       `1.602176634e-19` (*probed*).
     - **Never put a hex pattern in a `const`.** It is an integer, and `CONST * x` then does
       integer arithmetic with only a warning (*probed*: the result was
       `0x4000000000000000`). A full-width hex pattern in a typed `f64` slot keeps its bits
       with no warning (*probed*), but the tables use decimals for review.
     - ±∞ are written `(1.0 / 0.0)`. NaN can't be a compile-time f64, so a NaN in a table
       is a generator error.
   - **Layout.** Data modules are **baked typed column arrays**, one per field
     (`var TM_<X>_<FIELD>: f64|u8|bool[N] = { .. };`), one element per line with a label
     comment.
     - Baked arrays are in the image, need no initialization, and are immune to
       `alloc_reset` (6.6.16).
     - Arrays of structs are refused, so a table becomes parallel columns.
   - **Verified.**
     - The generated files pass `cyrius fmt --check` and `cyrius lint` as written.
     - `scripts/oracle/probe_tables.py` compiles every constant and table cell through the
       Cyrius compiler and compares it bit for bit with the golden rows: 129 values,
       0 failures, 0 warnings. A tampered value fails.
     - `--check` covers the generated `.cyr` files as well as the corpus.
   - **The remaining data modules are generated at the start of their milestone.** Each
     one's column layout is part of that module's Cyrius design: `nucleus_data` and
     `decay_data` in M2, the timekeeping tables in M5, `atomic_data` in M3, `reaction_data`
     in M4.
     - Most tables have a public accessor, so they can be enumerated through the API: the
       lookups by (Z, A) and by Z, and the `Vec`-returning lists.
     - Three have none: the proton/neutron level orders, the leap-second table and the USNO
       segments. Nor do the model coefficients (SEMF, Myers–Swiatecki, AME2020's u→keV).
       **Decided in M2 (2026-11-08):** they are extracted mechanically by
       `scripts/oracle/src/rust_src.rs` (literal `const` items only, with checked lengths) and
       then **proven**: the emitter recomputes the oracle's output from each value over the whole
       grid and requires identical bits. A 1-ulp change to one coefficient fails the proof. M5's
       leap-second table and USNO segments follow the same path. Nothing is re-typed.
   - Class-A golden rows (§4) pin every generated constant and table cell bit for bit.
7. ✅ **Parity harness. Done 2026-11-08.**
   - **`tests/parity.tcyr`** replays every module listed in `tests/golden/modules.txt`. The
     generator writes that list, so a new module can't go untested. It checks each row
     within its class (README *Comparison*) and gives it one outcome:
     - ok;
     - **FAILED**: the first 10 per module are printed with the row's note and both bit
       patterns;
     - **measured**: class C with no bound set yet, recorded but not gated until step 8;
       outliers beyond 16 ulp or with a NaN/∞ mismatch are printed;
     - **pending**: the function isn't ported yet.
   - **Report.** Per function it prints the row counts and an ulp histogram (0, 1, 2, 3–16,
     >16, NaN/∞), which is the step 8 baseline.
   - **Exit code.** It is the number of failing modules, following the ecosystem repro
     convention, at most 13. A missing golden file fails its module.
   - **Dispatchers.** One per module maps a golden `fn` name to the Cyrius call. A function
     gets its entry the moment it is ported. So far:
     - `constants`: all 34, bit-identical;
     - `libm`: sqrt, fabs, floor, trunc and cbrt are class A, and all are bit-identical;
       sin, cos, tan, exp, log, sinh, pow and hypot are class C, measured. `libm::round`
       (the half-away helper, M5) and `libm::ldexp` (the exp2 scaling, M2) are pending.
   - **Buffers.** The token and statistics buffers are heap-allocated: as statics they
     would cross the compiler's 128 KiB static-data advisory, and a warning fails the gate.
   - **`tests/ulp.cyr`** (`lib/assert.cyr` has no float assertion) provides:
     - `ulp_bits_match`, `ulp_distance` (ordered-integer, ±0 as one point), `ulp_rel_match`
       and `ulp_rel_err`;
     - the assert wrappers `assert_bits_eq`, `assert_ulp` and `assert_rel`, which count
       through `assert`, so the reference-value tests use the same file.
   - **Verified.**
     - Negative tests all fail as they should, with exit 1: a tampered constant, a
       malformed token and a missing golden file.
     - `cyrius test`, `test --bounds` and `audit` pass, with no warnings.
   - **First x86_64 result:** 4,574 rows ok, 0 failed, 6,044 measured, 7,835 pending.
     - cbrt (2,074 rows), sqrt, sin, cos and tan are bit-identical to libm.
     - exp, log, sinh and pow are within 1 ulp, and hypot within 2.
     - **One real defect:** ganita `pow(-1, y)` gives +∞ for y > 2^64, and +0 for
       y < −2^64, instead of 1. It is filed as
       `ganita/docs/development/issues/2026-11-08-pow-minus-one-huge-even-exponent.md`, with
       a self-proving repro. tanmatra's own `pow` bases are all positive.
8. ✅ **Baseline measurements. Done 2026-11-08: [`port-parity.md`](port-parity.md).**
   - **Targets.** `libm.tsv` was run through every Cyrius counterpart on x86_64 ELF (native),
     x86_64 PE (wine) and aarch64 (qemu user mode). macOS isn't runnable here and is recorded as
     unmeasured.
   - **Bit-identical to Rust libm on every target:** sqrt, fabs, floor, trunc, **cbrt (all
     2,074 rows)**, sin, cos and tan. **exp** is bit-identical on aarch64; on x86_64 it is
     1 ulp off on 60 of 710 rows (x87). log, sinh and pow are ≤ 1 ulp, and hypot ≤ 2 ulp.
     Windows PE matches x86_64 ELF row for row.
   - **Class-C bound for the libm primitives: K = 2 ulp** on every measured target. Each side
     documents < 1 ulp, so two implementations can be 2 apart. Every observation is within it.
     A gating check with exp's K = 0 failed exactly x86_64's 60 rows.
   - **One accepted divergence:** `pow(-1, DBL_MAX)`, ganita's |x| == 1 defect. It is filed
     upstream and accepted by
     [ADR 0007](../adr/0007-accept-ganita-pow-minus-one-divergence.md).
   - **Conditioning finding.** `invariant_mass` is ill-conditioned in |p| when p ≫ m: a 1-ulp
     change in hypot moves it by more than 16 ulp on 103 of 515 rows. So **the port
     reproduces Rust libm's hypot as `_tm_hypot`**, and the four-momentum functions become
     class B. A prototype, given in port-parity.md, matched all 600 hypot rows bit for bit on
     x86_64 and aarch64.
   - **Rule for tanmatra's own class-C bounds.** K_F = 2 · S · (class-C calls on its path) + 2,
     where S is F's measured sensitivity, set before F's first comparison. If K_F > 64, make
     the primitive exact instead.

**Gate:**
- `cyrius build` and `cyrius bench` pass. `bench` proves the step-2 manifest fix.
- `cyrius test` runs the parity harness (0 modules yet) with no warnings.
- `cyrius lint` and `cyrius fmt --check` are clean.
- Both generators and `tests/ulp.cyr` are committed, along with the upstream issues.

---

## 3. Milestones

Each milestone ends with the same gates:
- `cyrius build`. `--strict` has been a no-op since 6.3.2, because undefined functions are
  already errors.
- `cyrius test` and `cyrius test --bounds`, run on the module's golden rows plus its ported
  reference-value tests, and **with no warnings**.
  - The f64 type diagnostics come only from the compiles that `test` and `build` drive, not
    from `cyrius lint`. Each one is a bug: *integer arithmetic with an f64 right operand*,
    *f64 arithmetic with a non-f64 right operand*, *an integer stored into an f64 slot*,
    *comparison mixes f64 and integer*.
  - From M6, also `cyrius test --features optics`.
- `cyrius bench`, with the module's benches recorded against the Rust numbers.
- `cyrius vet <file>`, `cyrius lint <file>…` and `cyrius fmt --check <file>…`, which take file
  arguments, and `cyrius audit`, which runs fmt, lint, test and bench over the project.
- `cyrius distlib --all`, then `cyrius distlib --check` for the base bundle and every profile;
  then `cyrius deps --locked`.
- **Surface parity** (ranga pattern):
  - every in-scope Rust `pub fn` in the module has a Cyrius counterpart or a written reason;
  - `rust-old/src/X.rs` `pub fn` names are normalized and diffed against `src/X*.cyr` `fn`
    names;
  - the same diff is the M7 release gate.
- CHANGELOG and `state.md` updated.

| Milestone | Modules | Key work | Golden grid (minimum) |
|---|---|---|---|
| **M1 Foundation** | error, constants, particle, relativity | Generated CODATA 2022 `const`s. `TanmatraError` → a code enum plus an `err_out` detail. `FourMomentum` as an inline **struct value** with four typed f64 fields and a `FourMomentum_add` overload for `+` (Rust implements only `Add`). `momentum_magnitude` on a bit-exact **`_tm_hypot`** (port-parity.md; it moves the `libm::hypot` dispatch and class A with it). Enum properties (`Quark`, `Lepton`, `Boson`, `FundamentalForce`) read the generated `particle_data.cyr` columns, with a bounds check on the enum value. Probe the `Nucleus` value struct inside a `Result` (§5) | every enum variant; β ∈ {0, ±0.5, ±(1−1e-12), ±1, ±1.5}; four-momenta from 1 keV to 10 GeV |
| **M2 Nuclear** | nucleus_data, nucleus, shell_model, decay_data, decay | Generated AME2020, radii, moments, superallowed and 114-isotope **baked column tables**. SEMF + Myers–Swiatecki (`ganita_f64_cbrt`). Mayer–Jensen and empirical p/n orders; Brennan–Bernstein coupling. `decay_chain` with the IT ground-state rule. **Bateman scaling-and-squaring matrix exponential**, with `ldexp` as a multiply by `f64_exp2(k)` | all (Z, A) with 1 ≤ Z ≤ 118, Z ≤ A ≤ 3Z (binding, shell-corrected, J^π); all isotopes; the U-238, U-235 and Th-232 chains; Bateman with equal λ, the 9-member long-lived chain, and the U-238 series at 1 y, 10³ y and 4.463e9 y |
| **M3 Atomic** | atomic_data, atomic_config, atomic_hydrogenic, atomic_qed | Generated IE (118) and EA (118, 3-state) tables; `ElectronAffinity` as a `: stack` enum. Configurations with the 20 NIST exceptions. Rydberg and reduced-mass lines. Dirac and first-order fine structure. Radial functions (Laguerre coefficients and factorial-ratio normalization) for n ≤ 60. Exact radial dipole integrals → Einstein A/B. Lamb shift with Bethe logarithms, VP, hyperfine, g-factors, Breit. Write the `powi` sites in Rust's evaluation order (§5) | Z = 1..118 (configs, IE, EA); n1 < n2 ≤ 12, Z ∈ {1, 2, 3, 26, 92}; R_nl for all n ≤ 10, l < n, at 50 radii; A-values for every allowed n_u ≤ 10 → n_l; all (n, l, j) with n ≤ 4 |
| **M4 Reactions & scattering** | reaction_data, reaction, scattering | ENDF/B-VIII.0 thermal cross sections, resonance integrals and yields (generated). `NuclearReaction` layout with an optional projectile, neutron counts and a lepton vec. Presets. Rutherford (`ganita_f64_tan`), Mott (p·β), Born, the form factor (series branch), Klein–Nishina (series branch), Maximon pair production (both branches). The 7-argument `breit_wigner_cross_section` and `born_screened_coulomb_with_masses` keep their Rust signatures | all presets; energies on both sides of each branch point (KN γ = 1e-2, form factor qR = 0.5, pair k = 4); θ ∈ {1e-3, 0.1, π/2, π} |
| **M5 Timekeeping** | timekeeping | CIPM 2025 standards. Leap-second table, plus the 1961–1971 USNO segments and MJD → calendar. `AtomicInstant` as a value struct in integer nanoseconds; `+|` at Rust's three `saturating_add` sites; local round-half-away (§5); `nanoseconds_since` returns a `(seconds, nanos)` tuple. Geoid-referenced satellite correction. Sagnac (loop and one-way). TCG/TT and TCB/TDB conversions | leap lookups 1961–2030; `AtomicInstant` add/diff at epoch 0, at 2.1e9 s and at the ±9.2e18 s edges, plus **exact .5 ns ties** in `from_tai_seconds` and `add_seconds`; J2000 and T0 conversions |
| **M6 Optics profile** | optics (prakash 2.15.1) | `[features] optics` + optional `[deps.prakash]` + `[lib.optics]`. `lines_to_spd` → `spd_new`; the Storey & Hummer Balmer series; A-value series intensities; `line_to_rgb` → `wavelength_to_rgb` | Balmer SPD at 0.1 and 1 nm FWHM; invalid-range guards |
| **M7 Release 2.0.0** | all | JSON for the exchange types (`Nucleus`, `Isotope`, `FourMomentum`, `AtomicInstant`, `NuclearReaction`) if a consumer needs them. Use `#derive(Serialize)` only once cyrius issue `2026-11-07-derive-serialize-key-scan-reads-past-nul` is fixed; until then, write it by hand on bayan. Rust-test ledger (mela pattern): each of the 476 + 20 + 2 Rust tests is ported, folded into golden rows, or marked not applicable / deferred-bridge. Full surface-parity diff. `cyrius distlib --all`. `docs/benchmarks-rust-v-cyrius.md`. `cyrius audit`. VERSION 2.0.0 | full golden corpus green on every CI target |

Suggested order: M1 → M2 → M5 → M3 → M4 → M6 → M7. Timekeeping is self-contained and
exercises the integer-time paths early. Atomic is the largest module.

---

## 4. Parity policy

Every golden row carries one tolerance class:

| Class | Tolerance | Applies to |
|---|---|---|
| **A — exact** | bit-identical, or equal integers/strings | constants and table cells; configurations; J^π; decay chains; leap seconds; `AtomicInstant`; magic numbers; enum properties |
| **B — arithmetic** | ≤ 2 ulp, **expected 0** | formulas using only + − × ÷ √, `ganita_f64_cbrt` and `_tm_hypot`: SEMF, Rutherford, Klein–Nishina (series branch), relativity including the four-momentum functions, Q-value arithmetic. Both sides round each operation correctly and neither contracts to FMA. A non-zero distance means the operation order differs from Rust, so check that first |
| **C — transcendental** | ulp distance ≤ K. For the libm primitives K = 2; for each tanmatra function, K is derived by the [port-parity rule](port-parity.md#rule-for-the-class-c-bounds-of-tanmatras-own-functions) before its first comparison | anything through exp/ln/pow/sin/cos/tan/sinh: shell term, Dirac, Lamb shift, Mott, pair production, form factor, radial functions |
| **D — iterative** | relative ≤ 1e-12 | Bateman matrix exponential, radial dipole integrals, numerically integrated references |

Rules:
- **Never widen a class to make a row pass.** A failure beyond its class is a defect in one of
  three places:
  - in the port;
  - in the Cyrius math function: file it upstream, as was done for `f64_cbrt`;
  - in the Rust oracle: fix Rust, bump to 1.3.x and regenerate.
- **The one exception is a divergence accepted by an ADR** in `docs/adr/`. An example is a row
  where Cyrius is correctly rounded and `libm` is not (ranga's ADR 002 is the precedent). The
  harness lists such rows by name, and they are still reported.
- **The physics reference tests travel with the code.** The audit's reference-value tests are
  ported as Cyrius assertions as well as golden parity: NIST A-values, measured Lamb shifts,
  AME2020 binding energies, IERS offsets and ENDF values. Parity proves "same as Rust"; the
  reference tests prove "right".
- **Branch points are tested on both sides:**
  - the Klein–Nishina series switch;
  - the form-factor series switch;
  - the pair-production k = 4 join;
  - the degenerate-λ case in Bateman;
  - the pre/post-1972 TAI−UTC switch;
  - the .5 ns rounding tie.

---

## 5. Rust → Cyrius mapping (tanmatra specifics)

### Math functions

This table covers every `libm` call in tanmatra 1.3.0 (`grep -oh 'libm::[a-z0-9_]*' src`) and
every core f64 method it uses.

| Rust (calls) | Cyrius | From | Parity |
|---|---|---|---|
| `sqrt` (23), `fabs` (9), `floor` (2), `trunc` (1) | `f64_sqrt`, `f64_abs`, `f64_floor`, `f64_trunc` | compiler builtins | exact |
| `sin` (15), `cos` (7) | `f64_sin`, `f64_cos` | builtins, needing `math` on every target | fdlibm, ≤ 1 ulp, same bits on every target; class C |
| `exp` (9), `log` (6) | `f64_exp`, `f64_ln` | builtins, needing `math` on aarch64 | ≤ 1 ulp; x87 on x86_64, so the bits are **not** the same across targets; class C per target |
| `cbrt` (7) | `ganita_f64_cbrt` | ganita 1.2.15 | correctly rounded, as is `libm` 0.2.16's CORE-MATH `cbrt` → bit-identical; class B |
| `pow` (5) | `ganita_f64_pow` | ganita (fdlibm `e_pow`) | ≤ 1 ulp, exact for representable integer powers; class C. Keep `x^(5/3)` as `x·cbrt(x)²`, as in Rust |
| `sinh` (1), `tan` (1) | `ganita_f64_sinh`, `ganita_f64_tan` | ganita | class C. `tan` no longer needs a local sin/cos (ganita 1.2.11) |
| `hypot` (2: `momentum_magnitude`) | private `_tm_hypot`: Rust libm's (musl) algorithm, operation for operation | tanmatra (M1) | **bit-identical** (prototype 600/600 on x86_64 and aarch64; [port-parity.md](port-parity.md)). Not `ganita_f64_hypot`, which is ≤ 2 ulp from libm, too far for the ill-conditioned four-momentum functions |
| `round` (2: `timekeeping.rs:523`, `:600`) | private `_tm_round_half_away`: Rust libm's `trunc(x + copysign(0.5 − 0.25ε, x))` (M5), bit-identical on all 659 `libm::round` rows and every .5 ns tie | tanmatra | ⚠ `f64_round` rounds **ties to even** (2.5 → 2, −2.5 → −2; *probed*). `libm::round` rounds ties away from zero. A precedent exists: abaco `src/dsp.cyr` `f64_round_half_away` |
| `ldexp` (2: Bateman scaling) | private `_tm_ldexp`: Rust libm 0.2.16's scalbn transcribed (prescale by 2^1023 or 2^-969, then 2^n from bits) | tanmatra (M2) | **bit-identical**: 306/306 `libm::ldexp` rows (class A), including the subnormal and overflow edges. There is no ldexp/scalbn in the stdlib |
| `f64::powi` (`atomic.rs:782`, `:832`, `:1107`) | explicit products in LLVM's square-and-multiply order: `x*x`, `x*(x*x)`, `(x*x)*(x*x)` | — | bit-identical to Rust. `((x*x)*x)*x` is **not** |
| `f64::midpoint` (2, for j + ½) | `(a + b) / 2.0` | — | identical for in-range arguments |
| `i128` (`AtomicInstant::nanoseconds_since`) | the i128 itself as a pair `(hi, lo)`: the signed high and low 64 bits, computed by 128-bit subtract with borrow, × 10⁹ through `mulh64`, and a sign-extended add with carry (M5). For results that fit an i64, hi is the sign and lo the value | — | exact at the ±9.2e18 s and i64 MIN/MAX edges (432 golden rows, class A). There is no i128, and `u128` supports only + and − and can't be a parameter, return or field type |
| `f64 as i64` / `as u32` casts | `_tm_f64_to_i64_sat` (M5): truncation, NaN → 0, saturation at the i64 range, as Rust's `as` | tanmatra | the TAI − UTC calendar arithmetic and `from_tai_seconds` |
| `saturating_add` (`timekeeping.rs:503`, `:580`, `:606`) | `+|` | stdlib `overflow` | `+?`, `-?` and `*?` **abort the process with exit 57** (*probed*), so they are forbidden in `src/` |

Some functions are absent from the 6.7.18 stdlib: fma, expm1, log1p, log10, copysign, fmod
and isnan. tanmatra uses none of them; test for NaN with `x != x`.

**Bit patterns.** A bare copy `var u = x;` of an f64 is itself f64-typed (since 6.6.11), so `u < v` is a *float*
compare. To work on the bits, take them through a helper with an untyped parameter and an `: i64` return:
`fn _tm_bits(x): i64 { return x; }`. The first hypot prototype missed `hypot(±inf, NaN)` on exactly this.

**Typing rule.** Bind every call result to a `: f64` local before using it in an operator.
- The ganita functions return an untyped `i64`.
- Binary operators are typed by their **left** operand.
- An untyped left operand silently does integer arithmetic. `var q = ganita_f64_pow(…); q * 2.0`
  gives 0 (*probed*).

### Names (M1)

- **Free functions** keep their Rust names (`lorentz_gamma`).
- **A Rust enum's methods** become functions of the enum value, prefixed with snake(enum):
  `Quark::mass_mev(q)` is `quark_mass_mev(q)`. Their values are read from the generated
  `*_data.cyr` tables.
  - Rust can't hold an out-of-range enum value, but a Cyrius caller can. It gets NaN from an f64
    property, −1 from an integer one and 0 from a predicate, and never an abort.
- **A Rust struct** becomes `struct Tm<Name>` with an `impl Tm<Name>` block, so the methods keep
  their names.
- **Globals and enum members** carry `TM_`, and private helpers `_tm_`.
- **Return types.** Enum and `u8` returns are refused, so integer and predicate results are
  `: i64` (0/1).
- No M1 name collides with the 6.7.18 stdlib or prakash (scanned).
- `scripts/surface_parity.py` checks every ported Rust item against these rules.

### Data model

| Rust | Cyrius |
|---|---|
| `Nucleus { z: u32, a: u32 }` (Copy, validated) | **A packed `u64` word `z \| (a << 32)`**: the struct's own bytes, with `nucleus_z(n)` / `nucleus_a(n)` / `nucleus_n(n)` accessors. `nucleus_new` validates and returns a value-form `Result`. **Decided by M1's probe (2026-11-08):** `Ok(n)` of an 8-byte struct carries its bytes, but binding the payload back to a struct (`var n: S = f()?`) reads it as an address and segfaults, silently. That is filed as `cyrius/docs/development/issues/2026-11-08-small-struct-result-payload-read-as-address.md`. A plain word rides `Result` and `?` safely |
| `FourMomentum` (4 × f64, Copy, `impl Add`) | **`struct TmFourMomentum`** with typed f64 fields and an `impl TmFourMomentum` block (M1). Constructors are `TmFourMomentum_new` / `_at_rest` / `_from_mass_and_momentum`, methods are called `p.invariant_mass()`, and `fn add(self, rhs)` is the `+` operator. Build it in a local and return the local, because `return S{…}` is refused |
| `AtomicInstant { seconds: i64, nanos: u32 }` | **`struct TmAtomicInstant { seconds: i64; nanos: i64; }`**, 16 bytes by value, with an `impl` (M5): `TmAtomicInstant_new` (nanos taken mod 2^32 as Rust's u32, then carried), `_from_tai_seconds`, and methods `tai_seconds`, `to_tt_seconds`, `to_gps_seconds`, `duration_since`, `nanoseconds_since`, `add_nanoseconds`, `add_seconds`. All arithmetic is integer, `+|` saturates the seconds, and the nanosecond carry is Euclidean |
| `#[non_exhaustive] enum` (DecayMode, OrbitalType, TimeScale, FrequencyStandard, …) | `enum` with explicit, append-only numeric values. Property functions use `match` with an `_ =>` arm, because `match` warns on a non-exhaustive enum. Case values may be any const expression |
| `ShellLevel { n_shell, l, two_j }` (Copy) | `struct TmShellLevel { n_shell; l; two_j; }` with `impl` `new`, `degeneracy`, `j` and `label` (a `Str`), returned by value. `shell_model_levels()` is `shell_model_level_count()` + `shell_model_level(k)`. `shell_occupation` fills a caller `TmShellLevel[32]` and an occupation buffer and returns the count. **Changed after M3:** M2 packed it into a word `n_shell \| (l << 21) \| (two_j << 42)`, but Rust's public u32 fields allow a hand-built level past 2²¹, where the word silently printed the wrong label and degeneracy. Golden probes up to u32::MAX now pin it |
| `NuclearMoments` from `Option<NuclearMoments>` | `struct TmNuclearMoments` written through an out-pointer: `nucleus_nuclear_moments(n, &out)` returns 1/0. A 16-byte struct can't ride a value-form `Option` (M2) |
| `SuperallowedDecay` | `struct TmSuperallowedDecay { parent; daughter; ft_seconds: f64 }` with `impl` `corrected_ft_seconds() → Option`. `superallowed_ft_values()` is `superallowed_ft_value_count()` + `superallowed_ft_value(k)` (M2) |
| `Option<f64>` / `Option<u32>` returns | the value-form `Option` (`lib/tagged.cyr`): `var t, v = f(..); is_some(t)`, with the f64 bound by `var x: f64 = v;` (M2) |
| `(u32, i32)` / `(f64, f64)` returns | a multi-value `return (a, b);`, read as `var a, b = f(..);` (M2: `ground_state_spin_parity`, `superallowed_average_ft`) |
| `ElectronAffinity::{Bound(f64), Unbound, Unknown}` | **Decided in M3:** `enum TmElectronAffinity { BOUND; UNBOUND; UNKNOWN; }`, the kind. `electron_affinity(z)` returns `Result` `Ok(kind)` and the bound value is `electron_affinity_ev(z)` (0.0 otherwise, as in Rust). The planned `enum ElectronAffinity: stack { Bound(ev); .. }` can't be used: `electron_affinity` returns `Result<ElectronAffinity>`, and a `: stack` enum is itself one (tag, payload) pair, so it can't nest inside a value-form `Result` |
| `QuantumNumbers { n, l, ml, ms }` from `QuantumNumbers::new → Result` | `struct TmQuantumNumbers` with `impl` `new(n, l, ml, ms, out) → Result` `Ok(out)`, writing through the out-pointer (M3). Four u32/i32 fields don't fit a packed word, and a 32-byte struct can't ride a `Result` |
| `OrbitalFilling { n, orbital, electrons }` and `electron_configuration → Vec` | `struct TmOrbitalFilling { n; orbital; electrons; }` in a caller array: `var cfg: TmOrbitalFilling[20]; electron_configuration(z, &cfg, 20)` returns `Ok(count)`, and entry k is `base + k` with `var base: *TmOrbitalFilling = &cfg;` (a typed pointer steps by `sizeof`). ⚠ `&cfg + k` is an untyped address and steps by **bytes** (*probed*). `format_configuration(cfg, len)` and `format_configuration_short(cfg, len, z)` return a `Str`. **Decided in M3:** a 21-bit packed word, as `ShellLevel` is, failed the golden rows for a hand-built entry with n = electrons = 2³⁰, which Rust's public u32 fields allow. The struct round-trips every u32 |
| The named series (`lyman_series` … `pfund_series`, `Vec<(u32, f64)>`) | caller buffers: `lyman_series(z, n_max, n_out, nm_out, cap)` returns `Ok(count)` and writes up to `cap` (n, λ) pairs; the count needs no loop past `cap` (M3) |
| Integer arguments outside Rust's u32 / i32 | not a value Rust can be called with: a fallible function returns its `Err`, an infallible f64 one NaN, a `TransitionType` one `FORBIDDEN` (M3, `atomic_config.cyr`). Where Rust's u32 arithmetic would overflow (`2n − 1` for n ≥ 2³¹: a panic under the oracle's overflow checks), the port computes exactly in i64 |
| `Result<T, TanmatraError>` | value-form `Result` with **`?`** propagation (*probed*). The Err payload is a `TM_ERR_*` code (`error.cyr`). **Decided in M2 (2026-11-08): no detail travels with the code.** Every Rust error string only restates the caller's own input (the half-life, Z and A), and `tm_error_message(code)` gives the fixed text, so there is no `err_out`. Never return a `Result` through `callptr`: its payload is lost |
| Identity structs carrying strings or vecs (`Isotope`, `NuclearReaction`, `SuperallowedDecay`, `ThermalCrossSection`, `FissionYield`) | heap structs with typed fields, read as `p.field` through a `*T`. **Avoid `#derive(accessors)` getters on f64 fields.** They are untyped, so `P_x(p) + P_y(p)` is an integer add |
| `Vec<T>` returns (configurations, shell occupation, series, yields, chains, Bateman populations) | `lib/vec.cyr` of pointers, or f64 buffers with an explicit length. Hot APIs take a caller buffer |
| `known_isotopes()` (allocates 114 structs per call) and every data table | baked typed column arrays from M0 step 6, with indexed accessors. There is no lazy initialization and no `alloc_reset` hazard. M2: `known_isotope_count()` + `known_isotope(k)` returns a `struct TmIsotope` by value, whose name is a C string in a baked `u8` table with offsets. That is faster than Rust's `Vec` (0.82×) |
| `decay_chain(start, max)` / `bateman_chain(λ, N₀, t)` (Vec returns) | caller buffers with a count: `decay_chain(start, max, nuclei_out, modes_out, cap)` and `bateman_chain(λ_ptr, n, N₀, t, out)`. Bateman's n×n scratch is `fl_alloc`'d and `fl_free`'d inside, so a consumer's loop doesn't grow the heap; -1 if the scratch can't be allocated. A non-finite λt, where Rust's scaling loop never terminates, stops at 2200 squarings (M2) |
| constants | generated decimal `const TM_<NAME>`s (M0 step 6). **No hex in a `const`** |
| `String` names and labels | `Str` via `str_builder` |

### Cyrius rules that shape tanmatra's code

- **Names.** A bundle is concatenated into its consumers, so every public global carries the
  `TM_` prefix:
  - constants, as `TM_ELECTRON_MASS_MEV`;
  - enum members, as `TM_QUARK_UP`, inside `enum TmQuark`;
  - table columns, as `TM_QUARK_MASS_MEV`.

  The golden token `e:Quark.WPlus` maps to `TM_QUARK_`… by the rule `TM_` + snake(enum) +
  `_` + snake(variant) (`scripts/oracle/README.md`). Private helpers use `_tm_`.
- **Arguments.** Arguments past the sixth are passed on the stack, so the two 7-argument
  functions keep their signatures. Named and default arguments (6.7.7) keep call sites
  readable.
- **Negation.**
  - Write `-x` or `0.0 - x`.
  - `(0 - x)` on an f64 is an integer subtraction of the bit pattern (*probed*; it warns).
  - Negative literals such as `-1.5` and `-0.0` work.
- **Literals.** There is no exponent notation. The generator expands every value to positional
  decimal, or emits hex inside a typed `f64` slot.
- **Arrays and scratch memory.**
  - Array locals have been per-call stack storage since 6.3.15. Arrays above ~120 KB fall back
    to static storage, with a warning.
  - Never return an array local.
  - Bateman and radial-integral scratch uses a caller buffer or an arena. The heap is still a
    bump allocator, so a long-running consumer loop must not leak through it.
- **Globals.** The 4,096 non-literal-initializer cap was removed in 6.6.9.
- **Unsigned types.**
  - `u8`, `u16` and `u32` exist as locals, parameters (narrowed) and packed fields, but not as
    return types.
  - `u64` behaves as `i64`: there is no unsigned compare.
  - Validate non-negative Z, A, n, l and j at every API boundary.
- **Shifts.** `>>` is logical and `>>>` is arithmetic, the reverse of Java and JavaScript.
  Choose by the signedness of the Rust operand.
- **Struct values** changed heavily between 6.7.10 and 6.7.17.
  - Keep the shapes simple: no struct-valued calls at top level, and pass a `*T` when the
    callee mutates.
  - Re-probe struct values at every pin bump.
- **SIMD.** Don't use the `f64v_*` builtins until cyrius issue
  `2026-11-06-f64v-builtins-no-scalar-tail` is fixed: they overrun `dst` when `n` isn't a
  multiple of the lane count.

---

## 6. After 2.0.0

1. **Bridges** (`src/bridge.cyr`), ported per consumer as each lands in Cyrius:
   - chrono first, against stdlib `lib/chrono.cyr`, which is unblocked now;
   - kimiya, bijli and bhava;
   - jyotish, hisab-mimamsa and falak;
   - kiran/joshua (`SimulationClock`).

   Parity comes from the Rust bridge tests in `rust-old/`.
2. **soorat integration** (`src/soorat.cyr`), once soorat ports.
3. **Removing `rust-old/`.** This happens only after three things:
   - every deferred bridge and integration is ported or explicitly dropped;
   - the golden corpus and the generators are archived;
   - the M7 Rust-test ledger is closed.

   prakash removed its `rust-old/` in 2.2.3 and mela in 1.1.6; both keep it recoverable from git.
4. **Roadmap items continue in Cyrius:**
   - the v1.1 remainder: Y_lm, oscillator strengths, Hartree screening, He⁺/Li²⁺ validation;
   - v1.6: simulation support.

---

## 7. Tracking

| Milestone | Status | Golden rows | Reference tests | Benches | Notes |
|---|---|---|---|---|---|
| M0 Initialize | ✅ done 2026-11-08 | 18,453; 10,617 ok, 0 failed, 1 accepted (ADR 0007), the rest pending until ported | — | — | Steps 0–8 done; the generators grow per milestone by design. Parity is green on x86_64 ELF, x86_64 PE and aarch64, with libm class C bounded at 2 ulp ([port-parity.md](port-parity.md)). Filed upstream: the `cyrius port` scaffold issue and the ganita pow issue. Next: M1, which starts by adding `_tm_hypot` and the dispatchers |
| M1 Foundation | ✅ done 2026-11-08 | particle 102 + relativity 6,768 + constants 34: **all bit-identical (0 ulp)** on x86_64 ELF, x86_64 PE and aarch64 | 45 Rust tests ported (24 particle, 21 relativity; 5 serde ones don't port), plus 1 Cyrius-only: 114 assertions | 3 (`lorentz_gamma_1000` 2.4× Rust; `invariant_mass_1000` 6.8×; `particle/properties_all` 19×) | `numerics.cyr` (`_tm_hypot`), `error.cyr`, `particle.cyr`, `relativity.cyr`; `[lib]` and `dist/tanmatra.cyr`; surface parity 76/76. Probe: `Nucleus` becomes a packed word (cyrius issue filed) |
| M2 Nuclear | ✅ **done 2026-11-08** (M2a nucleus and shell model; M2b decay and Bateman) | nucleus 74,312 (74,297 at M2, plus 15 large-field `ShellLevel` probes after M3) + decay 8,511: **all class A/B rows 0 ulp**; `remaining_fraction` (C, K = 4) ≤ 1 ulp; Bateman (D) ≤ 9 ulp on x86_64, **0 on aarch64**; `libm::ldexp` 306/306 bit-identical | nucleus.tcyr (58 tests + 1 assert), decay.tcyr (see state.md) | 8 (nucleus 4, decay 4: see benchmarks doc) | `nucleus_data.cyr` and `decay_data.cyr` generated, with the private values proven; `nucleus.cyr`, `shell_model.cyr`, `decay.cyr`; `_tm_ldexp`; surface parity 46/46 for M2; the error-detail design decided (none) |
| M5 Timekeeping | ✅ done 2026-11-08 | timekeeping 20,703: **all 0 ulp** on all three targets (leap lookups 1955–2035, 14,695 TAI − UTC rows, `AtomicInstant` at the ±9.2e18 s and MIN/MAX edges with the 128-bit `nanoseconds_since`, exact .5 ns ties); `sagnac_correction_ns` (C, K = 6) 0 ulp; `libm::round` 659/659 bit-identical. **The corpus has 0 pending rows** | timekeeping.tcyr: 69 Rust tests (193 assertions); 5 not carried over (4 serde, 1 deprecated) | 1 (`atomic_instant_add_seconds_1000` 2.6×) | `timekeeping_data.cyr` generated: the leap table derived from the API and proven, the USNO segments extracted and proven; `timekeeping.cyr`; `_tm_round_half_away`, `_tm_f64_to_i64_sat`; surface parity 43/43 |
| M3 Atomic | ✅ done 2026-11-08 | atomic 43,309: **every class A/B row 0 ulp** on all three targets, and Einstein A/B (D) **0 ulp** too. Lamb shift (C, K = 10) ≤ 2 ulp on x86_64, ≤ 1 on aarch64; R_nl (C, K = 6) ≤ 2 and R²r² (C, K = 18) ≤ 7 on x86_64, **0 on aarch64**. Dirac is class B (sqrt only), not C | atomic.tcyr: 110 Rust tests (597 assertions); the 4 serde ones are not carried over | 7 (1.9–3.3× Rust for the computed bodies, 6–43× for the call-bound ones; see the benchmarks doc) | `atomic_data.cyr` generated: IE and EA through the API with each entry's citation; the filling order and Bethe logarithms extracted, the 20 NIST exceptions derived from `electron_configuration`, all proven. `atomic_config.cyr`, `atomic_hydrogenic.cyr`, `atomic_qed.cyr`; surface parity 41/41; `validate.py` re-derives the configurations and the class-B lines independently |
| M4 Reactions & scattering | ☐ | | | | |
| M6 Optics profile | ☐ | | | | |
| M7 Release 2.0.0 | ☐ | | | | |
| Bridges (deferred) | ⏸ | | | | chrono is unblocked (stdlib `lib/chrono.cyr`). The rest wait on bijli, kimiya, bhava, jyotish, falak, hisab-mimamsa and kiran/joshua |
| soorat integration (deferred) | ⏸ | | | | waits on soorat |
