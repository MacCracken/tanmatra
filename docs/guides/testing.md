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
| Unit tests | 191 | `src/*.rs` (`#[cfg(test)]` modules) |
| Integration tests | 20 | `tests/integration.rs` |
| Doc tests | 1 | `src/lib.rs` |
| **Total** | **212+** | |

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

Current coverage: **91.7%** (895/976 lines).

## Benchmarks

11 criterion benchmarks in `benches/benchmarks.rs`:

```bash
cargo bench
```

| Benchmark | Typical Time |
|-----------|-------------|
| `nucleus/binding_energy_1000` | ~18 us |
| `atomic/spectral_line_1000` | ~1.6 us |
| `atomic/electron_config_36` | ~2.0 us |
| `decay/decay_chain_10` | ~3.5 us |
| `nucleus/shell_occupation_126` | ~13 us |
| `atomic/ionization_energy_118` | ~66 ns |
| `relativity/lorentz_gamma_1000` | ~4.6 us |
| `scattering/rutherford_1000` | ~7.4 us |
| `decay/bateman_chain_3` | ~88 ns |
| `atomic/radial_wavefunction_100` | ~1.5 us |
| `decay/known_isotopes_alloc` | ~3.0 us |
