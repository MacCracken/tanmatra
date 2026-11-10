# Development Roadmap

**Forward-facing only.** This lists what is not done yet. Shipped work is in [`CHANGELOG.md`](../../CHANGELOG.md),
the port's milestones (M0–M7, 2.0.0 and 2.0.1) are in the [port plan](tanmatra-port-plan.md)'s tracking table, and
live status is in [`state.md`](state.md). A finished item is deleted from here, not checked off.

## How releases are picked

- **2.0.x, patches:** fixes, upstream follow-ups, new tests of existing functions, and work that just missed a
  release.
- **2.x, minors:** new public functions. Each keeps the whole 2.0 API, so a consumer upgrades without changes.
- **3.0.0:** the only release that may break the 2.x API. No breaking change is planned. The backlog below holds the
  research-scale items, which have no design yet.

Work that ports Rust code has its parity reference in `rust-old/` and gets golden rows from `scripts/oracle/`, as the
port did. New physics has no Rust counterpart, so there are no golden rows to match. Instead, each new function is
pinned by reference-value tests against a cited primary source (CODATA, NIST, PDG, NNDC). A fix that changes what Rust
1.3.0 returned needs an ADR first.

## 2.0.x — patches

None of these adds a function.

| Item | Waits on | Change |
|---|---|---|
| **He⁺ and Li²⁺ lines (Z > 1 validation)** | nothing | Reference-value tests of `spectral_line_vacuum_nm(z, a, n1, n2)` at ⁴He⁺ and ⁷Li²⁺, against the NIST ASD He II and Li III lines. The bound is set beforehand from what the formula leaves out (fine structure, QED). A defect the tests find is fixed here, with an ADR, since the fix diverges from Rust |
| **ganita `pow(−1, y)` for \|y\| > 2⁶⁴** | `ganita/docs/development/issues/2026-11-08-pow-minus-one-huge-even-exponent.md` | Once the cyrius pin moves to a stdlib with the fix, the one accepted libm row comes back ok. Then remove it from `_par_accepted` in `tests/parity.tcyr` and mark [ADR 0007](../adr/0007-accept-ganita-pow-minus-one-divergence.md) *Superseded* |
| **chrono dates before 1970** | `cyrius/docs/development/issues/2026-11-09-chrono-dates-before-1970-silently-wrong.md` | Once `lib/chrono.cyr` is fixed, extend `tests/bridge.tcyr`'s chrono recipe back to 1958 and drop the "from 1970 on" caveat in `src/bridge.cyr`'s header. The bridge functions don't change |

A cyrius or prakash pin bump ships as a patch once the full cleanliness check passes on x86_64 ELF, x86_64 PE and
aarch64.

## 2.1.0 — atomic structure

The Rust roadmap's v1.1 remainder, plus the two bridge helpers that wait on no consumer. Nothing blocks any of it.

- **Spherical harmonics Y_lm**, the angular part of the hydrogenic wavefunction. The radial part R_nl is already
  ported, to n = 60.
- **Oscillator strengths from first principles:** absorption f-values from the exact hydrogenic radial dipole
  integrals that the Einstein A and B coefficients already use. Reference values come from the NIST ASD.
- **Hartree screening:** effective nuclear charge for multi-electron atoms. Its reference data is chosen and cited
  before the code is written.
- **The prakash-labelled bridge helpers,** `transitions_to_wavelengths` and `nuclear_spin_to_hyperfine_scale`, added
  to `src/bridge.cyr`. These are ports, so their 4 Rust tests and their golden rows come from `rust-old/`. prakash is
  already in Cyrius. The port plan held them back for the bridge phase, which 2.0.1 opened.

## 2.2.0 — simulation support

The Rust roadmap's v1.6. All of it is new physics, tested against reference values.

- **Monte Carlo decay simulation:** stochastic chains over `decay_chain`, tested against Bateman's expectation.
- **Particle transport:** energy loss, range and straggling, against the NIST stopping-power and range tables (PSTAR,
  ASTAR, ESTAR).
- **Detector response functions:** Gaussian smearing.
- **Phase-space generators** for multi-body decays, on `TmFourMomentum`.

Monte Carlo needs a seeded, reproducible generator, and tanmatra has none. The stdlib's `lib/random.cyr` is the
kernel CSPRNG, which can't be seeded, and hisab 3.4.1 exports no generator either. The release's first decision is which
generator to use, recorded in an ADR, without adding a dependency just for it.

## Later 2.x — as consumers land in Cyrius

Each bridge waits for its consumer to be in Cyrius (port plan §6). When the consumer gets there, the bridge takes the
next free minor. Like the chrono bridge in 2.0.1, it lands in `src/bridge.cyr` with a golden grid from
`scripts/oracle/` and its Rust tests ported. Until then, `scripts/surface_parity.py` lists its items under
`DEFERRED_ITEMS`.

| Consumer | tanmatra surface | Rust tests | Blocked by |
|---|---|---|---|
| kimiya | `atomic_number_to_valence`, `neutron_count`, `mass_to_binding_deficit_mev`, `half_life_to_decay_constant` | 8 | kimiya 2.0.0. Its Cyrius port is under way (M3 of M7 done 2026-11-09) |
| bijli | `energy_to_wavelength_nm`, `nuclear_charge_field` | 4 | bijli has no repository yet |
| hisab-mimamsa | `gravitational_time_dilation`, `gravitational_time_offset_s` | 3 | not in Cyrius |
| falak | `tai_seconds_to_jd_tt`, `jd_tt_to_tai_seconds`, `tai_seconds_to_mjd_tt` | 4 | not in Cyrius |
| kiran / joshua | `SimulationClock` (7 functions) | 5, 1 of them serde | not in Cyrius |
| bhava | `TimeContext` (7 functions), built from a `SimulationClock` | 6, 2 of them serde | bhava is Rust (2.0.0). It also needs `SimulationClock` first |
| soorat | the `soorat-compat` integration as `src/soorat.cyr`, behind a feature as optics is: orbital density grids, nucleon positions, spectral lines and decay-chain graphs | 7 | soorat is Rust (1.0.0) |

**JSON for the exchange types** waits on a Cyrius consumer that exchanges them. The types are `Nucleus`, `Isotope`,
`FourMomentum`, `AtomicInstant` and `NuclearReaction`, with 7 Rust serde tests.
- Use `#derive(Serialize)` only once `cyrius/docs/development/issues/2026-11-07-derive-serialize-key-scan-reads-past-nul.md`
  is fixed. Until then, write the JSON by hand on bayan.
- A bridge type's serde tests port only if that type becomes an exchange type.

## Retiring `rust-old/`

This changes no API, so it ships in whichever release meets the three conditions (port plan §6):
- every deferred bridge and the soorat integration is ported or explicitly dropped;
- the golden corpus and its generators are archived;
- the [Rust-test ledger](rust-test-ledger.md) is closed. Today 48 tests are deferred: the 34 bridge tests, the 7 soorat
  tests and the 7 serde tests above.

If a consumer never reaches Cyrius, dropping its bridge is the way out. Dropping one diverges from Rust, so it takes
an ADR.

## 3.0.0 backlog

These items are research-scale and unscheduled. None has a design or a primary source yet, and each needs a scope
decision before it is pinned to a release. An item that fits the 2.x API moves up to a 2.x minor.

- Nuclear density functional theory interface
- Feynman diagram computation helpers
- Lattice QCD observable interfaces
- Real-time physics rendering with kiran, beyond the `SimulationClock` bridge
