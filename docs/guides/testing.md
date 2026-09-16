# Testing Guide

## Running Tests

```bash
# All tests with all features
cargo test --all-features

# Default features only
cargo test

# Specific module
cargo test --lib nucleus

# Specific test
cargo test --lib fe56_binding_energy

# With output
cargo test --lib -- --nocapture
```

## Test Categories

| Category | Count | Location |
|----------|-------|----------|
| Unit tests | 476 (all features), 464 (no default features) | `src/*.rs` (`#[cfg(test)]` modules) |
| Integration tests | 20 | `tests/integration.rs` |
| Doc tests | 2 | `src/lib.rs`, `src/timekeeping.rs` |
| **Total** | **498** | |

### Reference-Value Tests

Physics tests assert values from primary data with tolerances set by the
reference uncertainty, never ranges wide enough to pass a wrong formula:
AME2020 masses and binding energies, NUBASE2020 half-lives and J^π, NIST ASD
ionization energies and A-values, ENDF/B-VIII.0 cross sections and yields,
CODATA 2022 constants, measured Lamb shifts and hyperfine splittings, IERS
Conventions time-scale offsets, and exact results (normalization integrals,
atom conservation, numeric integration of differential cross sections).

## Test Patterns

### Serde Roundtrip

Every public type has a serde roundtrip test:

```rust
#[test]
fn serde_roundtrip_nucleus() {
    let n = Nucleus::iron_56();
    let json = serde_json::to_string(&n).unwrap();
    let back: Nucleus = serde_json::from_str(&json).unwrap();
    assert_eq!(n, back);
}
```

### Physics Validation

Tests compare computed values against known experimental/reference data:

```rust
#[test]
fn fe56_binding_energy_per_nucleon() {
    let fe56 = Nucleus::iron_56();
    let bea = fe56.binding_energy_per_nucleon();
    // Experimental: ~8.790 MeV; semi-empirical within ~2%
    assert!(bea > 8.6 && bea < 9.0);
}
```

### Error Path Coverage

All error conditions are tested:

```rust
#[test]
fn invalid_quantum_numbers() {
    assert!(QuantumNumbers::new(0, 0, 0, 1).is_err()); // n=0
    assert!(QuantumNumbers::new(1, 1, 0, 1).is_err()); // l >= n
}
```

### Exhaustive Enum Variant Coverage

All variants of each enum are exercised in at least one test to ensure coverage of match arms.

## Coverage

Target: **90%+**

```bash
# Using tarpaulin
cargo tarpaulin --all-features --skip-clean

# Using llvm-cov (if available)
cargo llvm-cov --all-features --html --output-dir coverage/
```

Coverage was last measured at 91.7% before the 2026-09-16 repairs; re-run
tarpaulin to refresh.

## Benchmarks

19 criterion benchmarks in `benches/benchmarks.rs` (timings from
`cargo bench -- --warm-up-time 1 --measurement-time 2`, 2026-09-16):

```bash
cargo bench
```

| Benchmark | Typical Time |
|-----------|-------------|
| `nucleus/binding_energy_1000` | ~16.8 µs |
| `nucleus/binding_energy_shell_corrected_1000` | ~129 µs |
| `nucleus/shell_occupation_126` | ~18.8 µs |
| `nucleus/ground_state_spin_parity_100` | ~0.84 µs |
| `atomic/spectral_line_1000` | ~1.04 µs |
| `atomic/electron_config_36` | ~1.93 µs |
| `atomic/ionization_energy_118` | ~61 ns |
| `atomic/radial_wavefunction_100` (n ≤ 3) | ~3.75 µs |
| `atomic/radial_wavefunction_n10_100` | ~8.2 µs |
| `atomic/einstein_a_to_n2_upto_n10` | ~7.8 µs |
| `atomic/lamb_shift_nlj_all_n4` | ~100 ns |
| `decay/decay_chain_10` | ~3.08 µs |
| `decay/bateman_chain_3` | ~1.0 µs (stable matrix exponential; was 88 ns with the unstable sum) |
| `decay/bateman_u238_series_15` | ~43 µs |
| `decay/known_isotopes_alloc` | ~2.32 µs |
| `relativity/lorentz_gamma_1000` | ~4.30 µs |
| `scattering/rutherford_1000` | ~6.77 µs |
| `scattering/pair_production_1000` | ~15.8 µs |
| `timekeeping/atomic_instant_add_seconds_1000` | ~11.6 µs |
