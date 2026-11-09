# scripts/oracle: what the port derives from the Rust oracle

`oracle` runs the frozen Rust tanmatra 1.3.0 in [`rust-old/`](../../rust-old/) and writes
two kinds of output:
- **The golden corpus, `tests/golden/<module>.tsv`** (port plan §2 step 5). It has one row
  per call, with every float as its IEEE-754 bit pattern. The Cyrius parity harness
  (`tests/parity.tcyr`, step 7) replays every row against the port and checks it within its
  tolerance class.
- **The generated Cyrius data modules, `src/constants.cyr` and `src/<module>_data.cyr`**
  (step 6). These are the constants and tables the port reads, so no physics value is ever
  typed by hand.

```sh
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml             # write both
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check  # verify; exit 1 on any difference
python3 -I scripts/oracle/validate.py .                                            # independent re-check of the corpus
python3 -I scripts/oracle/probe_tables.py .                                        # generated .cyr through the Cyrius compiler
```

## Rules

- **`rust-old/` is never built in place or edited.** This crate depends on it by path.
  `Cargo.lock` was seeded from `rust-old/Cargo.lock`, so every shared dependency is the
  version the oracle was tested with. `libm` is pinned `=0.2.16`.
- **The output is deterministic.** Running it twice gives byte-identical files, and
  `--check` passes on any machine. Grids use the pinned `libm` crate and seeded xorshift
  only, never the `std` f64 methods (they call the platform C library), timestamps or
  hash-map order.
- **The corpus is regenerated only from Rust 1.3.0.** A wrong oracle value is fixed in
  Rust: bump to 1.3.x, then regenerate. It is never patched in a `.tsv`.
- **Modules arrive with their milestone.** Each port milestone starts by adding its golden
  emitter (`src/m_<module>.rs`, with the grid from the port plan §3) and its data-module
  emitter (`src/t_<module>.rs`), and registering them in `src/main.rs`.
- **`--check` flags any `.tsv` in `tests/golden/` without an emitter as stale.**
- **`tests/golden/modules.txt`**, also generated, lists every module. `tests/parity.tcyr`
  replays exactly that list, so an emitter added here is replayed at once, as pending until
  its functions get dispatcher entries.

| Module | Milestone | Golden rows | Golden emitter | Generated Cyrius |
|---|---|---|---|---|
| `constants` | M1 | 34 | `src/m_constants.rs` | `src/constants.cyr` (`src/t_constants.rs`) |
| `particle` | M1 | 102 | `src/m_particle.rs` | `src/particle_data.cyr` (`src/t_particle.rs`) |
| `relativity` | M1 | 6,768 | `src/m_relativity.rs` | — (no tables) |
| `libm` | M0 step 8 | 11,549 | `src/m_libm.rs` | — |
| `nucleus` (with the shell model) | M2 | 74,312 | `src/m_nucleus.rs` | `src/nucleus_data.cyr` (`src/t_nucleus.rs`, with the private-value proofs) |
| `decay` | M2 | 8,511 | `src/m_decay.rs` | `src/decay_data.cyr` (`src/t_decay.rs`) |
| `timekeeping` | M5 | 20,703 | `src/m_timekeeping.rs` | `src/timekeeping_data.cyr` (`src/t_timekeeping.rs`: the leap table derived from `leap_seconds_at`, the USNO segments extracted, both proven) |
| `atomic` | M3 | 43,309 | `src/m_atomic.rs` | `src/atomic_data.cyr` (`src/t_atomic.rs`: the filling order and Bethe logarithms extracted, the 20 NIST exceptions derived from `electron_configuration`, all proven) |
| `reaction` | M4 | 7,965 | `src/m_reaction.rs` | `src/reaction_data.cyr` (`src/t_reaction.rs`: the presets, thermal cross sections, yields and pathways through the API; the resonance integrals and scattering's ζ(3) extracted and proven) |
| `scattering` | M4 | 10,704 | `src/m_scattering.rs` | — (its ζ(3) is in `reaction_data.cyr`) |
| optics | M6 | — | — | — |

## Format 1

### Header

Each file starts with five `#` lines:
- the format and module;
- the oracle provenance (tanmatra and libm versions);
- the regenerate command;
- the column list;
- the row count.

Lines starting with `#` are comments.

### Rows

```text
fn <TAB> class <TAB> inputs <TAB> outputs <TAB> note
```

- **`fn`** is the Rust path of what was called:
  - `Type::method` for methods and associated functions; the receiver is the first input;
  - the bare name for free functions;
  - `const` for a constant, whose input is its name;
  - `libm::<name>` for a math-library reference value.
- **`class`** is `A`, `B`, `C` or `D` (see *Comparison*).
- **`inputs` and `outputs`** are space-separated tokens, or `-` when empty.
  - A struct is flattened into its fields in declaration order. `FourMomentum` is four `f`
    tokens: energy, px, py, pz.
  - A tuple is flattened in order.
- **`note`** is the same call in readable form, e.g. `lorentz_gamma(0.5) -> 1.1547005383792517`.
  It is for review only: the harness ignores it.

### Tokens

| Token | Meaning |
|---|---|
| `f:<16 hex>` | f64, by IEEE-754 bit pattern (lowercase hex) |
| `i:<decimal>` | integer (u8 / u32 / i64 / usize); booleans are `0` / `1` |
| `e:<Enum>.<Variant>` | enum variant by its Rust name. In Cyrius, enum values follow the Rust declaration order and are append-only |
| `s:<text>` | UTF-8 string. Bytes outside printable ASCII (0x21–0x7e), and `%` itself, are `%XX` (uppercase hex) |
| `n:` | `Option::None`. `Some(x)` is just `x` |
| `E:<Variant>` | `Err(TanmatraError::<Variant>)`. Only the variant is compared, not the payload |
| `v:<len>` | a sequence of `len` elements, whose tokens follow in order. The harness knows each function's element shape |

A function returning a `Vec` too long for one row (`known_isotopes`) gets a row with its length (`v:<len>`) plus one row
per element under `fn[]` (`known_isotopes[]`), whose input is the index `i:k`.

### Comparison

**Non-float tokens** (`i`, `e`, `s`, `n`, `E`, `v` lengths) must be equal in every class.

**Floats** are compared by class:

| Class | A float output passes when |
|---|---|
| **A** exact | the bits are identical. Any NaN matches any NaN, because NaN payloads and signs differ between implementations |
| **B** arithmetic | the ulp distance is ≤ 2, measured on the ordered-integer line, where +0 and −0 are the same point. **Expected 0**: a non-zero distance means the operation order differs from Rust |
| **C** transcendental | the ulp distance is ≤ K, a bound per function (`_par_c_bound` in `tests/parity.tcyr`). K is derived and recorded in `docs/development/port-parity.md` before the function's first comparison. A function with no K yet is *measured*: its histogram is reported but not gated |
| **D** iterative | \|got − want\| ≤ 1e-12 · \|want\| |

In classes B, C and D:
- NaN matches only NaN, and ±∞ match only the same infinity;
- a zero `want` needs a zero `got`, of either sign (class D; in B and C, +0 and −0 are one point
  on the ulp line).

A tolerance is never widened to make a row pass (port plan §4). A divergence accepted by an
ADR in `docs/adr/` is listed by its exact function and input field in `_par_accepted`. It is
reported as `accepted (ADR n)` on every run.

## Generated Cyrius sources

Values come from the oracle's public API wherever one exists. The Rust source is also read
for its `///` doc comments, which carry the data citations (CODATA 2022, PDG 2024, ...) into
the Cyrius file.

### Private values (from M2)

Some values have no public accessor: the SEMF and Myers–Swiatecki coefficients, AME2020's
u→keV factor, `MS_MAGIC`, and the empirical proton and neutron level orders. They are never
re-typed. `src/rust_src.rs` extracts them mechanically under strict rules:
- it reads only a `const NAME: TYPE = <literal>;` item: a plain decimal, or a list of integer
  tuples (or `(integer, Enum::Variant)` pairs). The one other shape is a pure lookup function,
  `Some(match (a, b) { (1, 0) => <decimal>, .. _ => return None, })`, every line of which must be
  such an arm;
- a const's end is its first `;` outside a `//` comment (citations contain `;`);
- comments are read for labels only: each IE and EA entry's citation, checked to name its `(Z=k)`;
- a table's element count must match its declared length;
- floats are parsed with `str::parse::<f64>`, correctly rounded like rustc's own literal
  conversion.

**Every extracted value is then proven before anything is written.** The emitter transcribes
the Rust function that uses the value and recomputes the oracle's output from it, over the
whole milestone grid, requiring identical bits.
- `t_nucleus.rs` recomputes `binding_energy` and `binding_energy_shell_corrected` (all 14,160
  nuclei), `ground_state_spin_parity` (all 14,160), and every `experimental_atomic_mass_amu`.
- Nudging a single extracted coefficient by 1 ulp makes the generator refuse to write
  (*extraction proof failed: binding_energy(1, 2)*).

**M5.** The leap-second table needed no extraction: it is derived from the public `leap_seconds_at` (one entry per
change between consecutive months), then proven over years 1900–2200 and months 0–13. The 1961–1971 USNO segments
and `MJD_1972` are extracted, then proven by recomputing `tai_minus_utc_seconds_mjd` over 24,000 quarter-day MJDs
plus every segment start ±1e-6.

**M3.** The ionization energies and electron affinities come through the public API. `FILLING_ORDER` is extracted
and the 20 NIST exceptions are derived, as each element's two subshell changes from the pure Madelung filling;
together they are proven by recomputing `electron_configuration` for Z = 0–120 and u32::MAX. The Bethe logarithms
are extracted from `bethe_log`'s match and proven by recomputing `lamb_shift_nlj_ev` for Z = 0–118 and u32::MAX over
every (n ≤ 6, l ≤ n, 2j ≤ 2n + 1). A Bethe logarithm nudged by 1 ulp is refused (*Bethe-logarithm proof failed at
Z=1 n=3 l=0 2j=1*), and so is a dropped exception (*configuration proof failed at Z=45*).

**M4.** The presets, thermal cross sections, fission yields and nucleosynthesis pathways come through the public API
(their strings as NUL-terminated UTF-8 in one byte table). `RESONANCE_INTEGRALS` is extracted and proven by recomputing
all three lookups for every Z ≤ 130, A ≤ 300 and u32::MAX. Scattering's `ZETA3` is extracted and proven by recomputing
`pair_production_cross_section` on its k ≥ 4 branch for Z = 1–100 over 13 decades. `libm.tsv`'s sin, cos and log rows
are class A from M4: the port's `_tm_sin`, `_tm_cos` and `_tm_log` are Rust libm's algorithms.

### Emission rules

These were probed on Cyrius 6.7.18, and the generator enforces them.
- **f64s are decimal literals written positionally.** They use the shortest digits that
  round-trip, and each is checked to parse back to the oracle's bits; the lexer has no
  exponent notation. Decimal literals have been correctly rounded since 6.6.10.
- **Hex is never put in a `const`.** There it is an integer, and `CONST * x` does integer
  arithmetic.
- **±∞ are `(1.0 / 0.0)` and `-(1.0 / 0.0)`.** A NaN in a table is a generator error,
  because Cyrius refuses a compile-time NaN: its bits depend on the target.
- **Tables are baked column arrays.** Each is `var TM_<X>_<FIELD>: f64|u8|bool[N] = { .. };`,
  one element per line with a label comment. It lives in the image, needs no
  initialization, and is immune to `alloc_reset`.
- **The output is `cyrius fmt`- and `cyrius lint`-clean as generated.**

### Names

A bundle is concatenated into its consumers, so every generated global is prefixed:

| Oracle / golden token | Generated Cyrius |
|---|---|
| `const` row `s:ELECTRON_MASS_MEV` | `const TM_ELECTRON_MASS_MEV` |
| enum `Quark` | `enum TmQuark`, plus `const TM_QUARK_COUNT` |
| `e:Quark.Up`, `e:Boson.WPlus`, `e:TimeScale.TAI` | `TM_QUARK_UP`, `TM_BOSON_W_PLUS`, `TM_TIME_SCALE_TAI`: `TM_` + snake(enum) + `_` + snake(variant). snake() splits before a capital that follows a lowercase letter or digit, or that starts a new word after a run of capitals, so an acronym (`TAI`) stays one word. Values follow the Rust declaration order |
| `Quark::mass_mev` (tabulated) | `TM_QUARK_MASS_MEV[q]` |
| `FundamentalForce::mediator` (a slice) | `TM_FUNDAMENTAL_FORCE_MEDIATORS[START[f] .. START[f] + COUNT[f]]` |

### `probe_tables.py`

`probe_tables.py` compiles a program that prints every generated constant and table cell.
It builds in `build/oracle-probe/` with the project's manifest, then compares each value
bit for bit with the golden rows. It fails on any mismatch and on any compiler warning.

## `validate.py`

`validate.py` is an independent re-check of the corpus that doesn't use the Rust code. It
needs only Python 3, with no packages.
- **Format:** every row parses: five fields, a known class, well-formed tokens.
- **Constants:** `constants.tsv` matches the literals in `rust-old/src/constants.rs`.
- **Relativity:** the 2,962 class-A/B rows of the free functions and the arithmetic-only
  `FourMomentum` constructors and `add` are recomputed with Python floats in Rust's
  operation order, and must be bit-identical. Python floats are IEEE doubles.
- **libm:** `sqrt` matches the correctly rounded `math.sqrt`, and `cbrt` is exact on perfect
  cubes.

Extend it as modules arrive. A check that recomputes a row independently is worth more
than one that re-reads the oracle.
