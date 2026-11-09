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
| nucleus, shell model, decay | M2 | — | — | `nucleus_data.cyr`, `decay_data.cyr` |
| timekeeping | M5 | — | — | leap-second and frequency tables |
| atomic | M3 | — | — | `atomic_data.cyr` |
| reaction, scattering | M4 | — | — | `reaction_data.cyr` |
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

Values always come from the oracle's public API. The Rust source is read only for its `///`
doc comments, which carry the data citations (CODATA 2022, PDG 2024, ...) into the Cyrius
file. A table with no public accessor gets its emitter when its milestone designs the
Cyrius layout. That covers the proton/neutron level orders, the leap-second table and the
USNO segments: they are reached through the functions that use them, and they are never
re-typed.

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
| `e:Quark.Up`, `e:Boson.WPlus` | `TM_QUARK_UP`, `TM_BOSON_W_PLUS`: `TM_` + snake(enum) + `_` + snake(variant). Values follow the Rust declaration order |
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
