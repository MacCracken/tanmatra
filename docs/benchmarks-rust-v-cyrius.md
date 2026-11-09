# Benchmarks — Rust 1.3.0 vs the Cyrius port

Each row times one benchmark body: the same work in Rust (`rust-old/`, release build) and in
Cyrius (`tests/tanmatra.bcyr`, `cyrius bench`). Port milestones append rows as they land.
Numbers are host-specific, so compare within a column set, never across hosts.

## Method

- **Rust, criterion:** `cargo bench --manifest-path rust-old/Cargo.toml --bench benchmarks --
  <name>`, median of criterion's estimate. Where the Rust suite has no bench for a body, the
  same body is timed with `std::time::Instant` in a throwaway crate that uses `rust-old/` by path
  (median of 21 batches, `black_box` on every input and result). That harness gave
  `lorentz_gamma_1000` 4.63 µs against criterion's 4.78 µs, so the two methods agree to about
  3%.
- **Cyrius:** `cyrius bench`, average of 20 batches. lib/bench.cyr measures its clock-read floor
  and subtracts it from every sample (1.23 µs per read on this host).
- **Host:** x86_64 Linux, HPET clocksource. Toolchains: rustc 1.99.1 (release); Cyrius 6.7.18
  (ganita 1.2.15).

## M1 Foundation (2026-11-08)

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `relativity/lorentz_gamma_1000` | γ(i/1001), i = 1..1000 | 4.78 µs (criterion) | 11.52 µs | 2.4× |
| `relativity/invariant_mass_1000` | an electron at p = 1..1000 MeV/c along (3, 4, 12)/13: build the four-momentum, `invariant_mass()` | 21.2 µs (Instant) | 143.7 µs | 6.8× |
| `particle/properties_all` | 5 table-backed or computed properties for each of 6 variants | 15.1 ns (Instant) | 289 ns | 19× |

**Reading the gap.** Rust inlines these bodies into straight-line code: constant tables,
bounds checks gone, the hypot kernel inlined. Cyrius compiles each call as a call:
- every particle property makes a range check and a table read through a function call;
- `invariant_mass` makes two `_tm_hypot` calls, each with two `_tm_sq_lo` calls, plus the
  four-momentum by value.

The ratios grow with the share of call overhead in the body: lowest for `lorentz_gamma` (mostly
arithmetic), highest for the table lookups. Parity, not speed, is M1's bar, and every row above is
bit-identical to Rust. Candidates if a consumer needs the speed, each to be re-measured before
claiming it:
- `#inline` on `_tm_in_range` and `_tm_sq_lo`;
- the particle accessors reading the table directly after one range check.

## M2a Nuclear structure (2026-11-08)

All four are criterion benches in the Rust suite. The Cyrius bodies build the same nucleus lists outside the timed
loop.

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `nucleus/binding_energy_1000` | SEMF over Z = 1..100, A = Z..min(3Z, Z + 9) (984 nuclei) | 17.7 µs | 135.0 µs | 7.6× |
| `nucleus/binding_energy_shell_corrected_1000` | the same plus the Myers–Swiatecki term, A = Z..Z + 9 (1000) | 133.5 µs | 1.083 ms | 8.1× |
| `nucleus/ground_state_spin_parity_100` | odd Z = 1..99 at A = 2Z + 1 and 2Z | 0.955 µs | 13.6 µs | 14× |
| `nucleus/shell_occupation_126` | occupations for 1..126 nucleons | 13.5 µs (allocates a `Vec` per call) | 31.3 µs (caller buffers) | 2.3× |

**One optimization, measured.** `ground_state_spin_parity` first took 29.2 µs (30×): its level walk built a packed
level word through two calls for every level passed. It now reads only the 2j column and builds the word once, for
the level found, giving 13.6 µs. All 14,160 J^π golden rows are still exact.

**After M3, with `ShellLevel` a struct.** The J^π walk returns the level's index into the order tables and builds no
level value at all: 11.5 µs (12×). `shell_occupation` writes each level's three fields straight from the table instead
of packing a word through calls: 11.6 µs, against Rust's 13.5 µs (0.86×). The golden rows are unchanged.

## M2b Decay (2026-11-08)

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `decay/decay_chain_10` | the U-238 chain, at most 10 steps | 3.09 µs | 5.17 µs | 1.7× |
| `decay/bateman_chain_3` | λ = (0.1, 0.05, 0), N₀ = 10⁶, t = 100 s | 1.01 µs | 5.42 µs | 5.4× |
| `decay/bateman_u238_series_15` | the 15-member U-238 series at 1 y | 38.2 µs | 270.2 µs | 7.1× |
| `decay/known_isotopes_all` (Rust: `known_isotopes_alloc`) | all 114 isotope records | 2.36 µs (builds a `Vec` with a `String` each) | 1.93 µs (reads the baked table by value) | 0.82× |

The Cyrius U-238 series λ come from `decay_chain` and the isotope table, set up outside the timed loop; the Rust bench
types the same 14 half-lives inline.

**One optimization, measured.** Bateman's inner product first fetched both operands through a helper call per element:
574 µs (15×) for the U-238 series, 6.80 µs for `bateman_chain_3`. It now walks the two pointers directly, giving the same
products in the same order: 270 µs and 5.42 µs. All 28 Bateman golden rows give the same ulp histogram as before.

## M5 Timekeeping (2026-11-08)

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `timekeeping/atomic_instant_add_seconds_1000` | 1000 chained `add_seconds(1.000000001)` from (2.1e9 s, 123 456 789 ns) | 12.91 µs | 33.1 µs | 2.6× |

Each call does a trunc, a half-away round, two saturating integer adds and the Euclidean nanosecond carry, and
returns the 16-byte instant by value.


## M3 Atomic (2026-11-08)

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `atomic/spectral_line_1000` | every line n1 < n2 ≤ 10 for Z = 1..10 (450 calls) | 0.279 µs | 11.9 µs | 43× |
| `atomic/electron_config_36` | the configurations of Z = 1..36 | 1.77 µs (a new `Vec` each) | 3.32 µs (a caller array) | 1.9× |
| `atomic/ionization_energy_118` | every first ionization energy | 69.5 ns | 748 ns | 11× |
| `atomic/radial_wavefunction_100` | R₃₂(0.1·i a₀), i = 0..99 | 4.53 µs | 10.6 µs | 2.3× |
| `atomic/einstein_a_to_n2_upto_n10` | A for np → 2s and nd → 2p, n = 3..10 | 8.25 µs | 25.4 µs | 3.1× |
| `atomic/radial_wavefunction_n10_100` | R₁₀,₃(0.5·i a₀), i = 0..99 | 8.57 µs | 28.3 µs | 3.3× |
| `atomic/lamb_shift_nlj_all_n4` | the j = l + 1/2 shift of every (n, l) with n ≤ 4 | 103 ns | 640 ns | 6.2× |

**Reading the gap.** The bodies that compute (radial functions, Einstein A, configurations) are 2–3× Rust. The
call-bound ones are further off: `spectral_line_1000` is 450 calls of three divisions, which LLVM hoists and
pipelines across the loop (0.6 ns a call), while Cyrius makes each a call returning a `Result` pair (26 ns). The
ionization-energy lookup is the same story at 6 ns a call.

**Two optimizations, measured.**
- **`#inline` on the u32 / i32 range checks.** `_tm_u32_ok` and `_tm_i32_ok` are now branch-free single
  expressions, so `#inline` replays them at each call. `spectral_line_1000` went from 14.8 µs (53×) to 11.9 µs and
  `lamb_shift_nlj_all_n4` from 704 ns to 640 ns. Writing the checks out by hand reached 10.4 µs, which isn't worth
  the clutter.
- **The radial dipole integral, hoisted.** Rust recomputes pow(2/n, i), pow(2/n′, j) and the factorial product
  p!/c^(p+1) for every term of the double sum. The port computes the pows once per i and per j, and the products
  as one running product, so the calls and multiplications are the same and so are the bits. A Rust-shaped
  version timed 36.1 µs (4.4×) for `einstein_a_to_n2_upto_n10`; the hoisted one is 25.4 µs. All 2,886 Einstein
  golden rows stay 0 ulp.

## M4 Reactions & scattering (2026-11-08)

| Bench | Body | Rust | Cyrius | Cyrius / Rust |
|---|---|---|---|---|
| `scattering/rutherford_1000` | 5 MeV α on Au at θ = iπ/1001, i = 1..1000 | 6.55 µs | 62.2 µs | 9.5× |
| `scattering/pair_production_1000` | Pb at E_γ = 1.1 + i MeV, i = 0..999 (both Maximon branches) | 16.6 µs | 101.5 µs | 6.1× |

**Reading the gap.** Each Rutherford call is one sin plus a few operations, and sin dominates: Cyrius's `_tm_sin` takes
48 ns, while Rust's whole call, sin included, takes 6.5 ns. Pair production is mostly one log. The port uses its own
`_tm_sin` / `_tm_cos` / `_tm_log` because the port-parity rule requires bit-identical primitives for these formulas
(port-parity.md, M4). Measured on 1000 calls each, that exactness costs about 30%:

| Primitive | Cyrius builtin | The port's exact transcription |
|---|---|---|
| sin | `f64_sin`: 37.0 ns | `_tm_sin`: 48.1 ns |
| ln | `f64_ln` (x87): 37.8 ns; `_f64_ln_polyfill`: 56.4 ns | `_tm_log`: 53.8 ns |
