# Rust-test ledger

> **Generated** by `python3 -I scripts/rust_test_ledger.py .`; `--check` gates it (port plan M7, the mela
> pattern). Every test of the frozen Rust oracle (`rust-old/`, tanmatra 1.3.0) is ported to a Cyrius test,
> or recorded here with its reason. The golden corpus (`tests/golden/`) checks every ported function on its
> grid besides; this ledger is about the oracle's own assertions.

**498 Rust tests:** 424 ported, 1 partly ported, 25 not applicable, 48 deferred.

| Rust file | Cyrius tests | Rust tests | Ported | Partial | N/A | Deferred |
|---|---|---:|---:|---:|---:|---:|
| `src/particle.rs` | `tests/particle.tcyr` | 28 | 24 | 0 | 4 | 0 |
| `src/relativity.rs` | `tests/relativity.tcyr` | 22 | 21 | 0 | 0 | 1 |
| `src/error.rs` | — | 2 | 0 | 0 | 2 | 0 |
| `src/nucleus.rs` | `tests/nucleus.tcyr` | 63 | 58 | 1 | 3 | 1 |
| `src/decay.rs` | `tests/decay.tcyr` | 27 | 24 | 0 | 2 | 1 |
| `src/timekeeping.rs` | `tests/timekeeping.tcyr` | 74 | 69 | 0 | 3 | 2 |
| `src/atomic.rs` | `tests/atomic.tcyr` | 114 | 110 | 0 | 4 | 0 |
| `src/reaction.rs` | `tests/reaction.tcyr` | 43 | 39 | 0 | 3 | 1 |
| `src/scattering.rs` | `tests/scattering.tcyr` | 53 | 53 | 0 | 0 | 0 |
| `src/optics.rs` | `tests/optics.tcyr` | 5 | 4 | 0 | 1 | 0 |
| `src/bridge.rs` | `tests/bridge.tcyr` | 38 | 4 | 0 | 0 | 34 |
| `src/integration/soorat.rs` | — | 7 | 0 | 0 | 0 | 7 |
| `tests/integration.rs` | `tests/integration.tcyr` | 20 | 16 | 0 | 3 | 1 |
| doc tests (`src/lib.rs`, `src/timekeeping.rs`) | `tests/integration.tcyr` | 2 | 2 | 0 | 0 | 0 |
| **Total** | | **498** | **424** | **1** | **25** | **48** |

## Doc tests

- `rust-old/src/lib.rs`, the Quick Start example: `tests/integration.tcyr` `_t_doc_lib_quick_start`.
- `rust-old/src/timekeeping.rs`, AtomicInstant's example: `tests/integration.tcyr` `_t_doc_atomic_instant`.

## Tests without a full port

| Rust file | Test | Status |
|---|---|---|
| `src/particle.rs` | `serde_roundtrip_quark` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/particle.rs` | `serde_roundtrip_lepton` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/particle.rs` | `serde_roundtrip_boson` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/particle.rs` | `serde_roundtrip_force` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/relativity.rs` | `serde_roundtrip_four_momentum` | deferred: JSON for an exchange type, which M7 adds when a consumer needs it; no Cyrius consumer exists yet (port plan M7) |
| `src/error.rs` | `serde_roundtrip_error` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/error.rs` | `serde_roundtrip_error_with_string` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/nucleus.rs` | `serde_roundtrip` | deferred: JSON for an exchange type, which M7 adds when a consumer needs it; no Cyrius consumer exists yet (port plan M7) |
| `src/nucleus.rs` | `serde_roundtrip_shell_level` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/nucleus.rs` | `serde_roundtrip_nuclear_moments` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/nucleus.rs` | `corrected_ft_value_returns_average` | partial: its `corrected_ft_value` assert is n/a (deprecated in 1.3.0, not carried over, port plan §1); its `superallowed_average_ft` assert is ported as `_t_corrected_ft_value_returns_average_superallowed_part` |
| `src/nucleus.rs` | `serde_roundtrip_superallowed_decay` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/decay.rs` | `serde_roundtrip_decay_mode` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/decay.rs` | `serde_roundtrip_isotope` | deferred: JSON for an exchange type, which M7 adds when a consumer needs it; no Cyrius consumer exists yet (port plan M7) |
| `src/decay.rs` | `serde_roundtrip_isomeric_transition` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/timekeeping.rs` | `quality_factors_reasonable` | n/a: every assert calls `FrequencyStandard::quality_factor`, deprecated in 1.3.0 and not carried over (port plan §1); `quality_factor_for_linewidth` is ported and golden-checked |
| `src/timekeeping.rs` | `serde_roundtrip_frequency_standard` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/timekeeping.rs` | `serde_roundtrip_time_scale` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/timekeeping.rs` | `serde_roundtrip_atomic_instant` | deferred: JSON for an exchange type, which M7 adds when a consumer needs it; no Cyrius consumer exists yet (port plan M7) |
| `src/timekeeping.rs` | `atomic_instant_serde_normalizes` | deferred: JSON for an exchange type, which M7 adds when a consumer needs it; no Cyrius consumer exists yet (port plan M7) |
| `src/atomic.rs` | `serde_roundtrip_orbital_type` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/atomic.rs` | `serde_roundtrip_quantum_numbers` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/atomic.rs` | `serde_roundtrip_orbital_filling` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/atomic.rs` | `serde_roundtrip_transition_type` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/reaction.rs` | `serde_roundtrip_reaction` | deferred: JSON for an exchange type, which M7 adds when a consumer needs it; no Cyrius consumer exists yet (port plan M7) |
| `src/reaction.rs` | `serde_roundtrip_thermal_cross_section` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/reaction.rs` | `serde_roundtrip_nucleosynthesis` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/reaction.rs` | `serde_roundtrip_fission_yield` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/optics.rs` | `serde_roundtrip_spectral_line` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `tests/integration.rs` | `serde_roundtrip_nucleus` | deferred: JSON for an exchange type, which M7 adds when a consumer needs it; no Cyrius consumer exists yet (port plan M7) |
| `tests/integration.rs` | `serde_roundtrip_quark` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `tests/integration.rs` | `serde_roundtrip_decay_mode` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `tests/integration.rs` | `serde_roundtrip_quantum_numbers` | n/a: serde round-trip. JSON is carried only for the exchange types (port plan §1, M7) |
| `src/bridge.rs` | the other 34 | deferred: the bridges to bijli, kimiya, hisab-mimamsa, falak, kiran/joshua and bhava, and the prakash-labelled helpers, wait on their consumers being in Cyrius (port plan §1, *Deferred until the consumers exist*). The chrono bridge is ported |
| `src/integration/soorat.rs` | all 7 | deferred: the soorat integration waits on soorat (port plan §1) |
