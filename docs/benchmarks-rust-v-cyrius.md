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
