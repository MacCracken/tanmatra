# scripts/oracle: the golden corpus generator

`golden` runs the frozen Rust tanmatra 1.3.0 in [`rust-old/`](../../rust-old/) and writes
`tests/golden/<module>.tsv`: one row per call, with every float as its IEEE-754 bit
pattern. The Cyrius parity harness (`tests/parity.tcyr`, port plan §2 step 7) replays
every row against the port and checks it within its tolerance class.

```sh
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml             # write tests/golden/
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check  # verify; exit 1 on any difference
python3 -I scripts/oracle/validate.py .                                            # independent re-check (below)
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
- **Modules arrive with their milestone.** Each port milestone starts by adding its
  emitter (`src/m_<module>.rs`) with the golden grid from the port plan §3, and registering
  it in `modules()` in `src/main.rs`.
- **`--check` flags any `.tsv` without an emitter as stale.**

| Module | Milestone | Rows | Emitter |
|---|---|---|---|
| `constants` | M1 | 34 | `src/m_constants.rs` |
| `particle` | M1 | 102 | `src/m_particle.rs` |
| `relativity` | M1 | 6,768 | `src/m_relativity.rs` |
| `libm` | M0 step 8 | 11,549 | `src/m_libm.rs` |
| nucleus, shell model, decay | M2 | — | — |
| timekeeping | M5 | — | — |
| atomic | M3 | — | — |
| reaction, scattering | M4 | — | — |
| optics | M6 | — | — |

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
| **C** transcendental | \|got − want\| ≤ tol · \|want\|, where tol is per function and per target, measured in M0 step 8 from `libm.tsv` and recorded in `docs/development/port-parity.md` |
| **D** iterative | \|got − want\| ≤ 1e-12 · \|want\| |

In classes B, C and D:
- NaN matches only NaN, and ±∞ match only the same infinity;
- a zero `want` needs a zero `got`, of either sign.

A tolerance is never widened to make a row pass (port plan §4). A divergence accepted by
an ADR in `docs/adr/` is listed by row name in the harness.

## `validate.py`

`validate.py` is an independent re-check of the corpus that doesn't use the Rust code. It
needs only Python 3, with no packages.
- **Format:** every row parses: five fields, a known class, well-formed tokens.
- **Constants:** `constants.tsv` matches the literals in `rust-old/src/constants.rs`.
- **Relativity:** every class-A or class-B row is recomputed with Python floats in Rust's
  operation order, and must be bit-identical. Python floats are IEEE doubles.
- **libm:** `sqrt` matches the correctly rounded `math.sqrt`, and `cbrt` is exact on perfect
  cubes.

Extend it as modules arrive. A check that recomputes a row independently is worth more
than one that re-reads the oracle.
