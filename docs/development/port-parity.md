# Port parity — math-function baselines and the class-C bounds

Port plan §2 step 8: run the golden inputs through every Cyrius function tanmatra's math maps to, on every
target, against Rust libm 0.2.16. This page holds the ulp distributions, the documented guarantees on both
sides, the class-C bounds derived from them, and the rule for setting the bounds of tanmatra's own class-C
functions. The bounds live in `_par_c_bound` (`tests/parity.tcyr`), and this page is their rationale. Change
the two together.

**Measured 2026-11-08** with Cyrius 6.7.18 (ganita 1.2.15) against tanmatra 1.3.0's Rust libm 0.2.16, on
`tests/golden/libm.tsv`.

## Targets

| Target | How it ran | Measured |
|---|---|---|
| x86_64 Linux ELF | native (`cyrius test`) | ✅ |
| x86_64 Windows PE | wine (`cyrius test --win`) | ✅ identical to x86_64 ELF, row for row |
| aarch64 Linux ELF | qemu-aarch64 user mode (`cyrius test --aarch64`) | ✅ qemu emulates the FP instructions in IEEE softfloat; not yet confirmed on hardware |
| x86_64 / arm64 macOS | — | ❌ not runnable here. Expected to match x86_64 / aarch64 respectively (the same backends); confirm on a Mac before relying on it |

## Results

Per Cyrius counterpart: rows, then the ulp histogram of each output against Rust libm. Classes are from
`libm.tsv`.

| Rust libm | Cyrius | Class | Rows | x86_64 (ELF and PE) | aarch64 |
|---|---|---|---|---|---|
| `sqrt` | `f64_sqrt` | A | 774 | all 0 | all 0 |
| `fabs` | `f64_abs` | A | 374 | all 0 | all 0 |
| `floor` | `f64_floor` | A | 659 | all 0 | all 0 |
| `trunc` | `f64_trunc` | A | 659 | all 0 | all 0 |
| `cbrt` | `ganita_f64_cbrt` | A | 2,074 | all 0 | all 0 |
| `round` | port helper (M5) | A | 659 | pending | pending |
| `ldexp` | port helper (M2) | A | 306 | pending | pending |
| `sin` | `f64_sin` | C | 710 | all 0 | all 0 |
| `cos` | `f64_cos` | C | 710 | all 0 | all 0 |
| `tan` | `ganita_f64_tan` | C | 710 | all 0 | all 0 |
| `exp` | `f64_exp` | C | 710 | 650 × 0, **60 × 1** | all 0 |
| `log` | `f64_ln` | C | 874 | 864 × 0, 10 × 1 | 870 × 0, 4 × 1 |
| `sinh` | `ganita_f64_sinh` | C | 510 | 495 × 0, 15 × 1 | 509 × 0, 1 × 1 |
| `pow` | `ganita_f64_pow` | C | 1,220 | 1,206 × 0, 13 × 1, 1 accepted | same |
| `hypot` | `ganita_f64_hypot` | C | 600 | 503 × 0, 96 × 1, 1 × 2 | same |

Readings:

- **Correctly rounded on both sides means bit-identical.** sqrt is IEEE. For cbrt, Rust libm 0.2.16 has the
  CORE-MATH port and ganita has its own correctly rounded cbrt, and all 2,074 rows match. That includes every
  mass number 1–300 and every perfect cube up to 1000³. **Class B holds for cbrt, as the plan expected.**
- **Shared fdlibm kernels mean bit-identical in practice.** sin, cos and tan are fdlibm on both sides: 0 ulp
  on all 2,130 rows. pow is fdlibm `e_pow` on both sides, and 13 rows still differ by 1 ulp, in the
  extra-precision log2 path.
- **The targets differ only where the toolchain says they will.** On x86_64, `f64_exp` and `f64_ln` lower to
  x87 (`f2xm1`/`fyl2x`). On aarch64 they call `_f64_*_polyfill`, which is fdlibm-derived. So on aarch64,
  exp is bit-identical to Rust libm (which is fdlibm too), and log and sinh (built on exp) differ on fewer
  rows. Everything implemented in software in ganita (pow, hypot, tan, cbrt) is row-for-row identical across
  targets.

## Documented guarantees

| Function | Rust libm 0.2.16 (source comment) | Cyrius 6.7.18 / ganita 1.2.15 (source comment) |
|---|---|---|
| exp, log | FreeBSD msun: "the error is always less than 1 ulp" (`exp.rs`, `log.rs`) | `lib/math.cyr`: "finite results are within 1 ulp of the correctly rounded value", on every target, and "NOT promised bit-identical across targets" |
| sin, cos, tan | FreeBSD msun: "nearly rounded" (`tan.rs`) | `lib/math.cyr`: sin and cos within 1 ulp for every finite argument. ganita `tan`: fdlibm, "within 1 ulp everywhere it has been measured", worst found 0.940 ulp |
| pow | FreeBSD msun `e_pow`: "nearly rounded" | ganita: fdlibm `e_pow` core, "within 1 ulp everywhere", exact for representable integer powers |
| sinh | musl, built on expm1. No stated bound | ganita: built on expm1, which is measured at 0.87 ulp; correctly rounded where x³/6 is under half an ulp |
| hypot | musl (double-double squares, then sqrt). No stated bound | ganita: scaled `sqrt(1 + r²)`, which "doesn't overflow". No stated bound |
| cbrt | CORE-MATH: correctly rounded | ganita: "correctly rounded for every f64, subnormals included" |

## The class-C bounds

**libm primitives: K = 2 ulp**, on every measured target. Each side is within 1 ulp of the correctly rounded
value (or "nearly rounded", for fdlibm), so two of them can differ by up to 2 ulp. That is the bound the
guarantees support, and it is set from them, not from the sample.
- Every observation is within it: the maximum is 1 ulp, and hypot's is 2.
- hypot documents no bound on either side. Its measured maximum, 2 ulp, sits exactly at K, and a 3-ulp row
  would fail as a finding.
- The bounds are the same on all three targets, so `_par_c_bound` doesn't depend on the target. If a future
  target needs a different K, branch on it there and add a column above.

**One accepted divergence**, under [ADR 0007](../adr/0007-accept-ganita-pow-minus-one-divergence.md):
`libm::pow(-1.0, DBL_MAX)`, where ganita gives +inf and the answer is 1. This is ganita's missing
|x| == 1 guard, filed as `ganita/docs/development/issues/2026-11-08-pow-minus-one-huge-even-exponent.md`.
tanmatra never takes a negative base.

**The bounds gate.** With exp's K set to 0 as a check, x86_64 failed exactly its 60 one-ulp rows (exit 1),
and aarch64 passed all 710.

## Conditioning: why the four-momentum functions are class B, with an exact hypot

`FourMomentum::momentum_magnitude` is `hypot(hypot(px, py), pz)`. `invariant_mass`, `gamma`, `beta`,
`kinetic_energy` and `invariant_mass_two_body` sit on top of it. For p ≫ m, the mass is ill-conditioned in
|p|: δm/m ≈ (p/m)² · δp/p.

The 515 golden `invariant_mass` rows with finite, non-zero p and m were measured by moving |p| one ulp each
way and recording the change in the output:

| output change (ulp) | 0 | 1 | 2 | 3–16 | 17–10³ | 10³–10⁶ | > 10⁶ |
|---|---|---|---|---|---|---|---|
| rows | 337 | 41 | 6 | 28 | 36 | 30 | 37 |

103 rows move by more than 16 ulp. At the extreme, m flips between 0 and a positive value. ganita's hypot
differs from Rust's on 16% of rows, so on those rows **no class-C bound is meaningful**: the port would be
right, and every row would fail.

**Decision: the port reproduces Rust libm's hypot operation for operation, as a private `_tm_hypot`.** The
four-momentum functions then use only + − × ÷ and correctly rounded `sqrt`, so they are class B, expected
0 ulp. The golden generator now classes them B (`scripts/oracle/src/m_relativity.rs`); this tightens their
class rather than widening it.

A prototype of the M1 function matched **all 600 `libm.tsv` hypot rows bit for bit, on both x86_64 and
aarch64**. The prototype is the algorithm of Rust libm 0.2.16's `hypot.rs` (musl, MIT):

```
# The bit pattern of an f64 as an integer: an untyped parameter carries the word, and the
# `: i64` return types the caller's binding as an integer. A bare `var u = x;` copy of an
# f64 is f64-typed (6.6.11), so its `<` would be a FLOAT compare, false with a NaN: that
# made the first prototype miss hypot(±inf, NaN).
fn _tm_bits(x): i64 { return x; }

fn _tm_sq_lo(x: f64): f64 {            # the low part of x² (Dekker split at 2^27 + 1)
    var split: f64 = 134217729.0;
    var xc: f64 = x * split;
    var xh: f64 = (x - xc) + xc;
    var xl: f64 = x - xh;
    var hi: f64 = x * x;
    var lo: f64 = (((xh * xh) - hi) + ((2.0 * xh) * xl)) + (xl * xl);
    return lo;
}

fn _tm_hypot(xa: f64, ya: f64): f64 {
    var ux = _tm_bits(xa) & 0x7FFFFFFFFFFFFFFF;
    var uy = _tm_bits(ya) & 0x7FFFFFFFFFFFFFFF;
    if (ux < uy) { var t = ux; ux = uy; uy = t; }     # |x| >= |y|, by integer compare
    var ex = ux >> 52;
    var ey = uy >> 52;
    var x: f64 = ux;
    var y: f64 = uy;
    if (ey == 0x7FF) { return y; }                    # hypot(inf, nan) == inf
    if (ex == 0x7FF || uy == 0) { return x; }
    if (ex - ey > 64) { var s: f64 = x + y; return s; }
    var z: f64 = 1.0;
    var up: f64 = 0x6BB0000000000000;                 # 2^700
    var dn: f64 = 0x1430000000000000;                 # 2^-700
    if (ex > 0x3FF + 510) { z = up; x = x * dn; y = y * dn; }
    elif (ey < 0x3FF - 450) { z = dn; x = x * up; y = y * up; }
    var hx: f64 = x * x;
    var lx: f64 = _tm_sq_lo(x);
    var hy: f64 = y * y;
    var ly: f64 = _tm_sq_lo(y);
    var sum: f64 = ((ly + lx) + hy) + hx;
    var r: f64 = z * f64_sqrt(sum);
    return r;
}
```

When M1 lands `_tm_hypot`, two more steps follow:
- the harness's `libm::hypot` dispatch moves to it;
- `libm::hypot` becomes class A in `m_libm.rs`, which is again a tightening.

## Rule for the class-C bounds of tanmatra's own functions

Class C stays for formulas built on exp, ln, pow, sin, cos, tan or sinh, such as the shell term, Dirac, the
Lamb shift, Mott, pair production, the form factor and the radial functions. Each class-C function F gets
its bound **when it is ported, before its first comparison**:

1. **Measure F's sensitivity.** S is the largest output change, in ulp, produced by a 1-ulp change in any one
   class-C primitive result F uses, over F's golden grid.
   - It is found by perturbation, as above for `invariant_mass`. The script re-implements F's formula in
     Python, in Rust's operation order, with the primitive's result nudged with `nextafter`.
   - Alternatively, an analytic condition number may be used, with its derivation recorded.
2. **Set K_F = 2 · S · (class-C calls on F's longest path) + 2.** That is the primitives' 2 ulp carried
   through F, plus the final rounding.
3. **If K_F > 64, F is ill-conditioned on that grid.** Don't set a bound that large. Make the offending
   primitive bit-identical to the oracle's, as with hypot, which turns F into class B. If that isn't
   possible, write an ADR that splits the ill-conditioned rows out with their own documented bound.
4. **Record F, S, the method and K_F in the table below**, then set K_F in `_par_c_bound`.

Never raise K_F after seeing the port's results. A row that fails K_F is a defect, in the port, in Cyrius
(file it upstream), or in the Rust oracle (fix it, bump to 1.3.x, regenerate).

| Function | Milestone | S (method) | Class-C calls on its path | K_F |
|---|---|---|---|---|
| — none yet. The first are expected in M2 (the Myers–Swiatecki shell term) and M3 | | | | |

## Re-measuring

```sh
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check   # corpus current
cyrius test tests/parity.tcyr                 # x86_64
cyrius test --win tests/parity.tcyr           # x86_64 PE, needs wine
cyrius test --aarch64 tests/parity.tcyr       # aarch64, needs qemu-aarch64
```

Each run prints the per-function histograms that this page tabulates. Re-measure on every Cyrius pin bump.
The 6.7.x migration notes have changed math lowering before, and a changed histogram means a changed table
here: revisit the bounds deliberately, never by default.
