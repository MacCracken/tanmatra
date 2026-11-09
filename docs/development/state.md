# tanmatra — Current State

> Refreshed every release and at every port milestone. `CLAUDE.md` is
> preferences/process/procedures (durable); this file is **state** (volatile).

## Version

**2.0.0** (`VERSION`), the Cyrius library, 2026-11-08. The Rust crate's last release, 1.3.0, stays frozen
as the oracle in `rust-old/`. `cyrius.cyml` reads the version through `version = "${file:VERSION}"`, so
`dist/tanmatra.cyr` and `dist/tanmatra-optics.cyr` say v2.0.0. The tag and the release are the maintainer's.

## Port progress

Plan of record: [`tanmatra-port-plan.md`](tanmatra-port-plan.md).

| Milestone | Status |
|---|---|
| M0 Initialize | ✅ **done 2026-11-08.** Steps 0–4: Rust oracle frozen green; `cyrius port` scaffold; CLAUDE.md, roadmap and docs repaired; manifest configured; ganita 1.2.15 confirmed. **Step 5 built**: the golden generator `scripts/oracle/` writes `tests/golden/` (constants, particle, relativity, libm: 18,453 rows), deterministic, `--check`ed and independently validated. **Step 6 built**: the same tool generates `src/constants.cyr` (34 `TM_` consts) and `src/particle_data.cyr` (enums plus baked columns), all 129 values verified through the Cyrius compiler. **Step 7 done**: `tests/parity.tcyr` with `tests/ulp.cyr` replays the corpus: 4,574 rows ok, 0 failed, 6,044 measured (class C, no bound yet), 7,835 pending (not ported). **Step 8 done**: [`port-parity.md`](port-parity.md) has the libm baselines on x86_64 ELF, x86_64 PE and aarch64, with K = 2 ulp for the class-C libm primitives, one ADR-0007 accepted row, and the hypot conditioning finding (the port uses a bit-exact `_tm_hypot`). The harness now reports 10,617 ok, 0 failed, 0 measured, 7,835 pending, 1 accepted on all three targets |
| M1 Foundation | ✅ **done 2026-11-08.** `error`, `numerics` (`_tm_hypot`), `particle` and `relativity` are ported. Every particle, relativity and constants golden row is **0 ulp** on x86_64 ELF, x86_64 PE and aarch64. 45 Rust tests are ported (114 assertions), surface parity is 76/76, and 3 benches are recorded |
| M2 Nuclear | ✅ **done 2026-11-08.** M2a nucleus and shell model: all 74,312 golden rows 0 ulp (`ShellLevel` became a struct after M3, pinned by 15 large-field probes). M2b decay: all 8,511 rows pass. Class A/B are 0 ulp, `remaining_fraction` (C, K = 4) ≤ 1 ulp, and Bateman (D) ≤ 9 ulp on x86_64 and **0 on aarch64**. `libm::ldexp` is 306/306 bit-identical through `_tm_ldexp`. 82 Rust tests ported plus 1 assert (1,831 assertions); surface parity 46/46 for M2 |
| M5 Timekeeping | ✅ **done 2026-11-08.** All 20,703 timekeeping golden rows are 0 ulp on all three targets, including the 128-bit `nanoseconds_since` at the ±9.2e18 s and MIN/MAX edges and every .5 ns rounding tie. `libm::round` is 659/659, so **the corpus has 0 pending rows**. 69 Rust tests ported (193 assertions); surface parity 43/43 |
| M3 Atomic | ✅ **done 2026-11-08.** All 43,309 atomic golden rows pass on all three targets. Every class A/B row is 0 ulp (configurations, lines, fine structure, series, Zeeman/Stark, VP, Dirac, hyperfine, g-factors, Breit), and so are the 2,886 Einstein A/B rows (class D) on M3's grid (M6 found one transition, 49p → 1s, at ~2e-14 relative). The Lamb shift (C, K = 10) is ≤ 2 ulp on x86_64 and ≤ 1 on aarch64; R_nl (K = 6) ≤ 2 and R²r² (K = 18) ≤ 7 on x86_64, 0 on aarch64. 110 Rust tests ported (597 assertions); surface parity 41/41 |
| M4 Reactions & scattering | ✅ **done 2026-11-08.** All 18,669 reaction and scattering golden rows are 0 ulp on all three targets. The port-parity rule required bit-exact primitives (12 of 16 transcendental functions came out past K = 64), so `_tm_sin`, `_tm_cos` and `_tm_log` transcribe Rust libm; a 680,010-input sweep matches on every target. Every M4 function is class B except `rutherford_total_above_angle` (ganita tan, C, K = 10). 92 Rust tests ported (903 assertions); surface parity 26/26 + 18/18 |
| M6 Optics profile | ✅ **done 2026-11-08.** `[features] optics`, an optional prakash 2.15.1 and `[lib.optics]` (`dist/tanmatra-optics.cyr`); the core bundle is unchanged. All 1,039 optics rows pass on all three targets against the Cyrius prakash: `line_to_rgb` and the Balmer series 0 ulp, the SPDs ≤ 2 ulp (0 on aarch64), and 3 class-D series values at ~2e-14 relative (the 49p → 1s Einstein A). 4 Rust tests ported; surface parity 4/4. ⚠ The two optics test files carry prakash's three bundle warnings |
| M7 Release 2.0.0 | ✅ **done 2026-11-08.** The Rust-test ledger accounts for all 498 Rust tests: 420 ported, 1 partly, 25 n/a, 52 deferred (`docs/development/rust-test-ledger.md`, gated by `scripts/rust_test_ledger.py --check`). `tests/integration.tcyr` ports the integration and doc tests. Surface parity is the full diff. JSON for the exchange types is deferred (no Cyrius consumer yet). VERSION 2.0.0 |

## Toolchain

- **Cyrius pin**: `6.7.18` (in `cyrius.cyml` `[package].cyrius`)
- **ganita**: 1.2.15, from the 6.7.18 stdlib bundle (`lib/ganita.cyr`)
- **Rust oracle**: tanmatra 1.3.0 at `rust-old/`. It is green on rustc 1.99.1: 476 unit + 20 integration + 2 doc tests, clippy `-D warnings` and `cargo deny`

## Source

- **Rust reference:** 11,419 lines at `rust-old/`, frozen; do not edit. `rust-old/docs/` keeps the Rust testing guide and
  architecture overview; `rust-old/.github/` and `rust-old/scripts/` keep the Rust CI and the cargo scripts.
- **Cyrius port:** `[lib]` = `error` → `numerics` → `constants`* → `particle_data`* → `particle` → `relativity` →
  `nucleus_data`* → `nucleus` → `shell_model` → `decay_data`* → `decay` → `atomic_data`* → `atomic_config` →
  `atomic_hydrogenic` → `atomic_qed` → `reaction_data`* → `reaction` → `scattering` → `timekeeping_data`* →
  `timekeeping` (* generated by `scripts/oracle/`; never edit by hand). Bundled by `cyrius distlib` into `dist/tanmatra.cyr`. `src/main.cyr`
  is the scaffold stub.

## Tests

- `tests/tanmatra.tcyr`: scaffold smoke test (2 assertions)
- `tests/parity.tcyr`: the golden-corpus parity harness (below)
- `tests/particle.tcyr`: the 24 Rust particle tests, plus 1 Cyrius-only out-of-range test (80 assertions)
- `tests/relativity.tcyr`: the 21 Rust relativity tests (34 assertions)
- `tests/nucleus.tcyr`: 58 of the 63 Rust nucleus tests (1,765 assertions). Not carried over: the 4 serde tests, and
  the deprecated `corrected_ft_value` half of a fifth, whose `superallowed_average_ft` assert is ported
- `tests/decay.tcyr`: 24 of the 27 Rust decay tests (66 assertions); the 3 serde tests are not carried over
- `tests/timekeeping.tcyr`: 69 of the 74 Rust timekeeping tests (193 assertions). Not carried over: 4 serde tests, and
  `quality_factors_reasonable`, whose every assert calls the deprecated `quality_factor`
- `tests/tanmatra.bcyr`: benches; see `docs/benchmarks-rust-v-cyrius.md`
- `tests/atomic.tcyr`: 110 of the 114 Rust atomic tests (597 assertions); the 4 serde tests are not carried over
- `tests/reaction.tcyr`: 39 of the 43 Rust reaction tests (133 assertions); the 4 serde tests are not carried over
- `tests/scattering.tcyr`: all 53 Rust scattering tests (770 assertions)
- `tests/optics.tcyr`: the 4 Rust optics tests (8 assertions); the serde test is not carried over. Needs
  `cyrius deps --features optics`, like `tests/optics_parity.tcyr` and `tests/optics.bcyr`
- `tests/integration.tcyr`: the 16 non-serde tests of `rust-old/tests/integration.rs` and the 2 Rust doc tests (51
  assertions)
- `scripts/rust_test_ledger.py`: every Rust test ported or recorded with a reason (`docs/development/rust-test-ledger.md`);
  498 = 420 ported + 1 partial + 25 n/a + 52 deferred
- `tests/parity_core.cyr`: the replay machinery that `tests/parity.tcyr` (core modules) and
  `tests/optics_parity.tcyr` (optics) share
- `scripts/surface_parity.py`: every ported Rust public item has its Cyrius name. 254/254: error 6, constants 34,
  particle 18, relativity 18, nucleus 37, decay 9, timekeeping 43, atomic 41, reaction 26, scattering 18, optics 4
- `tests/tanmatra.fcyr`: scaffold fuzz stub
- Golden corpus, `tests/golden/*.tsv`:

  | Module | Rows |
  |---|---|
  | constants | 34 |
  | particle | 102 |
  | relativity | 6,768 |
  | nucleus | 74,312 |
  | decay | 8,511 |
  | timekeeping | 20,703 |
  | atomic | 43,309 |
  | reaction | 7,965 |
  | scattering | 10,704 |
  | optics | 1,039 |
  | libm | 11,549 |

  These are generated by `scripts/oracle/` (format spec in its README) and validated by `scripts/oracle/validate.py`,
  which also re-derives the atomic configurations and recomputes the class-B lines and levels independently.
- **`tests/parity.tcyr`**, the parity harness: 183,956 ok, 0 failed, 0 measured, 0 pending and 1 accepted (ADR 0007)
  on x86_64 ELF, x86_64 PE and aarch64; **`tests/optics_parity.tcyr`**, 1,039 ok on all three. Per-function histograms and the class-C bounds are in
  [`port-parity.md`](port-parity.md).

- **`tests/ulp.cyr`**: float comparison and assert wrappers for the parity harness and the reference tests.

Gates green on the scaffold (2026-11-08):
- `cyrius build`, `test`, `test --bounds` and `bench`
- `lint`, `vet` and `fmt --check`
- `deps --locked` and `audit`

## Dependencies

Direct (declared in `cyrius.cyml`):

- **stdlib:** syscalls, string, alloc, str, fmt, io, vec, math, ganita, overflow, tagged, result, freelist, sakshi
- **test stdlib** (`[test] stdlib`; test, bench and fuzz builds only): assert, bench
- **prakash 2.15.1:** optional, behind the `optics` feature (M6), commit-pinned in cyrius.lock with its hisab 3.4.1
  and sakshi 2.5.10. Its bundle needs the `sakshi` stdlib leaf, which `[deps]` lists

## Consumers

prakash, kiran and joshua consume the Rust crate. No Cyrius consumer builds against tanmatra yet.

## Upstream issues filed

- `prakash/docs/development/issues/2026-11-08-bundle-consumer-compile-warnings.md`. A consumer of the prakash 2.15.1 bundle can't compile warning-free:
  - two raw-bit compares warn (`lib/prakash.cyr:9213–9214`), with two checked silent fixes proposed;
  - hisab's `num_fft` is referenced but not shipped;
  - `sakshi` is both a stdlib leaf and a package.

  tanmatra **accepts** the three compile warnings, confined to `tests/optics.tcyr` and `tests/optics_parity.tcyr`,
  until prakash fixes them (decided 2026-11-08). The core suite is warning-free.

- `cyrius/docs/development/issues/2026-11-08-cyrius-port-clobbers-project-docs-and-omits-bench.md`. `cyrius port` overwrites
  CLAUDE.md, the roadmap, state.md and the workflows. It also pins `version = "0.1.0"` and omits `bench` from `[test] stdlib`.
  Its stubs fail `cyrius audit` (undocumented functions), and its `--dry-run` under-reports.

- `ganita/docs/development/issues/2026-11-08-pow-minus-one-huge-even-exponent.md`. `ganita_f64_pow(-1, y)` gives +∞
  (y > 2^64) or +0 (y < −2^64) instead of 1, with a self-proving repro (exit 12 on 1.2.15). tanmatra is unaffected: all its
  pow bases are positive.

- `cyrius/docs/development/issues/2026-11-08-small-struct-result-payload-read-as-address.md`. `Ok(s)` of a ≤ 8-byte struct
  carries its bytes, but `var s: S = <payload>` (or `= f()?`) reads it as an address and segfaults, silently. As a result,
  tanmatra's `Nucleus` is a packed word.

## Next

After 2.0.0 (port plan §6):
- the bridges, per consumer as each lands in Cyrius (chrono first, against stdlib `lib/chrono.cyr`);
- the soorat integration;
- JSON for the exchange types when a consumer needs it;
- the roadmap's physics items.

`rust-old/` stays until the deferred bridges and integrations are ported or dropped and the ledger is closed.
