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
  - `ElectronAffinity` becomes a **`: stack` payload enum**;
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
| `particle.rs` | 294 | 18 | 28 | `src/particle.cyr` |
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
   - `[lib]` order: error → constants → particle → relativity → nucleus_data → nucleus →
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
     - recomputes all 2,962 class-A/B relativity rows with Python floats in Rust's operation
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
6. **Table and constant generator in `scripts/oracle/`**, in the same crate as step 5. It
   reads the values through the oracle's public API, or from the primary source files the
   audit names, and never by parsing Rust source. Tables are generated, never re-typed: the
   audit found 17 of 32 AME rows and about 30 of 118 EA rows hand-transcribed wrong.
   - **Data modules** (`src/*_data.cyr`) are written as **baked typed column arrays**. Each
     array holds one field of the table, with a trailing comment giving the decimal value:
     ```
     var AME_Z: u16[32] = { 1, … };
     var AME_MASS_EXCESS_KEV: f64[32] = { 0x…, … };  # 7288.971064, …
     ```
     - Baked arrays are in the image, need no initialization, and are immune to
       `alloc_reset` (6.6.16).
     - Arrays of structs are refused, so a table becomes parallel columns.
     - Hex bit patterns are exact inside a typed `f64` slot (*probed*).
   - **`src/constants.cyr`** is generated from CODATA 2022 `allascii.txt` as decimal `const`s
     written out in full, e.g. `const ELEMENTARY_CHARGE_C = 0.0000000000000000001602176634;`.
     - The lexer has **no exponent notation**: `1.5e-3` is a parse error.
     - Decimal literals have been correctly rounded since 6.6.10. That literal gives
       `0x3c07a4da290c1653`, bit-identical to Python's `1.602176634e-19` (*probed*).
     - **Never put a hex pattern in a `const`.** It is an integer, and `CONST * x` then does
       integer arithmetic with only a warning (*probed*: the result was `0x4000000000000000`).
   - Class-A golden rows (§4) pin every generated constant and table cell bit for bit.
7. **Parity harness.**
   - `tests/parity.tcyr` loads `tests/golden/*.tsv` and checks each row within its class. It
     prints counts per module.
   - Its exit code is the number of failing modules, following the ecosystem repro
     convention. That is at most 13, so the 8-bit exit status can't wrap.
   - `lib/assert.cyr` has no float assertion, so add `tests/ulp.cyr` with
     `assert_bits_eq(got, want, label)` and `assert_ulp(got, want, max_ulp, label)`. The
     reference-value tests use the same file.
8. **Baseline measurements.** Before porting any formula, run the golden inputs through each
   Cyrius function tanmatra uses (§5 table) against Rust `libm`.
   - Run it **on every CI target**. On x86_64, `f64_exp` and `f64_ln` go through x87, so
     their bits differ from aarch64.
   - Record the ulp distributions per target in `docs/development/port-parity.md`. This sets
     class C (§4) from evidence, not guesswork.
   - Expected 0 ulp: `f64_sqrt` (IEEE) and `ganita_f64_cbrt` (correctly rounded on both
     sides). Any non-zero distance there is a defect.

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
| **M1 Foundation** | error, constants, particle, relativity | Generated CODATA 2022 `const`s. `TanmatraError` → a code enum plus an `err_out` detail. `FourMomentum` as an inline **struct value** with four typed f64 fields and a `FourMomentum_add` overload for `+` (Rust implements only `Add`). Enum methods (`Quark`, `Lepton`, `Boson`, `FundamentalForce`) as `match` tables. Probe the `Nucleus` value struct inside a `Result` (§5) | every enum variant; β ∈ {0, ±0.5, ±(1−1e-12), ±1, ±1.5}; four-momenta from 1 keV to 10 GeV |
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
| **B — arithmetic** | ≤ 2 ulp, **expected 0** | formulas using only + − × ÷ √ and `ganita_f64_cbrt`: SEMF, Rutherford, Klein–Nishina (series branch), relativity, Q-value arithmetic. Both sides round each operation correctly and neither contracts to FMA. A non-zero distance means the operation order differs from Rust, so check that first |
| **C — transcendental** | relative error ≤ the per-target bound measured in M0 step 8 (expected ~1e-14) | anything through exp/ln/pow/sin/cos/tan/sinh/hypot: shell term, Dirac, Lamb shift, Mott, pair production, form factor, radial functions |
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
| `sinh` (1), `tan` (1), `hypot` (2) | `ganita_f64_sinh`, `ganita_f64_tan`, `ganita_f64_hypot` | ganita | class C. `tan` no longer needs a local sin/cos (ganita 1.2.11) |
| `round` (2: `timekeeping.rs:523`, `:600`) | local `tm_round_half_away` | tanmatra | ⚠ `f64_round` rounds **ties to even** (2.5 → 2, −2.5 → −2; *probed*). `libm::round` rounds ties away from zero. A precedent exists: abaco `src/dsp.cyr` `f64_round_half_away` |
| `ldexp` (2: Bateman scaling) | multiply by `f64_exp2(k)` | builtin | exact for integer k in the normal range; there is no ldexp/scalbn in the stdlib |
| `f64::powi` (`atomic.rs:782`, `:832`, `:1107`) | explicit products in LLVM's square-and-multiply order: `x*x`, `x*(x*x)`, `(x*x)*(x*x)` | — | bit-identical to Rust. `((x*x)*x)*x` is **not** |
| `f64::midpoint` (2, for j + ½) | `(a + b) / 2.0` | — | identical for in-range arguments |
| `i128` (`AtomicInstant::nanoseconds_since`) | normalized `(seconds, nanos)` tuple | — | there is no i128, and `u128` supports only + and − and can't be a parameter, return or field type |
| `saturating_add` (`timekeeping.rs:503`, `:580`, `:606`) | `+|` | stdlib `overflow` | `+?`, `-?` and `*?` **abort the process with exit 57** (*probed*), so they are forbidden in `src/` |

Some functions are absent from the 6.7.18 stdlib: fma, expm1, log1p, log10, copysign, fmod
and isnan. tanmatra uses none of them; test for NaN with `x != x`.

**Typing rule.** Bind every call result to a `: f64` local before using it in an operator.
- The ganita functions return an untyped `i64`.
- Binary operators are typed by their **left** operand.
- An untyped left operand silently does integer arithmetic. `var q = ganita_f64_pow(…); q * 2.0`
  gives 0 (*probed*).

### Data model

| Rust | Cyrius |
|---|---|
| `Nucleus { z: u32, a: u32 }` (Copy, validated) | `struct Nucleus { z: u32; a: u32; }` **by value**: 8 bytes in one register, no allocation (*probed*). `nucleus_new` validates and returns a `Result`. If M1's probe finds that a struct payload won't ride a value-form `Result`, fall back to a packed `(z << 16) \| a` i64 with `nuc_z`/`nuc_a` accessors |
| `FourMomentum` (4 × f64, Copy, `impl Add`) | inline struct value with typed f64 fields and a `FourMomentum_add` overload for `+`. When returning by value, bind the literal to a local first, because `return S{…}` is refused |
| `AtomicInstant { seconds: i64, nanos: u32 }` | value struct `{ seconds: i64; nanos: u32; }` with a normalizing constructor |
| `#[non_exhaustive] enum` (DecayMode, OrbitalType, TimeScale, FrequencyStandard, …) | `enum` with explicit, append-only numeric values. Property functions use `match` with an `_ =>` arm, because `match` warns on a non-exhaustive enum. Case values may be any const expression |
| `ElectronAffinity::{Bound(f64), Unbound, Unknown}` | `enum ElectronAffinity: stack { Bound(ev); Unbound(); Unknown(); }`, which allocates nothing (*probed*). Bind the payload with `var e: f64 = ev;`. A user enum can't be a function's return type, so leave the function's return unannotated. A plain enum without `: stack` that carries a payload allocates |
| `Result<T, TanmatraError>` | value-form `Result` with **`?`** propagation (*probed*), plus a `TanmatraError` code enum. There is only one payload word, so detail strings go through an `err_out` pointer (the prakash pattern). Never return a `Result` through `callptr`: its payload is lost |
| Identity structs carrying strings or vecs (`Isotope`, `NuclearReaction`, `ShellLevel`, `SuperallowedDecay`, `ThermalCrossSection`, `FissionYield`) | heap structs with typed fields, read as `p.field` through a `*T`. **Avoid `#derive(accessors)` getters on f64 fields.** They are untyped, so `P_x(p) + P_y(p)` is an integer add |
| `Vec<T>` returns (configurations, shell occupation, series, yields, chains, Bateman populations) | `lib/vec.cyr` of pointers, or f64 buffers with an explicit length. Hot APIs take a caller buffer |
| `known_isotopes()` (allocates 114 structs per call) and every data table | baked typed column arrays from M0 step 6, with indexed accessors. There is no lazy initialization and no `alloc_reset` hazard |
| constants | generated decimal `const`s (M0 step 6). **No hex in a `const`** |
| `String` names and labels | `Str` via `str_builder` |

### Cyrius rules that shape tanmatra's code

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
| M0 Initialize | ◐ steps 0–4 done; step 5 built | 18,453 (M1 modules + `libm`) | — | — | 2026-11-08: the oracle is frozen green; the scaffold is in place and repaired; manifest and ganita 1.2.15 checked; all scaffold gates pass. The golden generator is in `scripts/oracle/`, deterministic and independently validated. Next: step 6 (tables and constants), step 7 (parity harness), step 8 (baselines from `libm.tsv`) |
| M1 Foundation | ☐ | | | | |
| M2 Nuclear | ☐ | | | | |
| M5 Timekeeping | ☐ | | | | |
| M3 Atomic | ☐ | | | | |
| M4 Reactions & scattering | ☐ | | | | |
| M6 Optics profile | ☐ | | | | |
| M7 Release 2.0.0 | ☐ | | | | |
| Bridges (deferred) | ⏸ | | | | chrono is unblocked (stdlib `lib/chrono.cyr`). The rest wait on bijli, kimiya, bhava, jyotish, falak, hisab-mimamsa and kiran/joshua |
| soorat integration (deferred) | ⏸ | | | | waits on soorat |
