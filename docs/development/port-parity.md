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
| `round` | **`_tm_round_half_away`** (M5): Rust libm's `trunc(x + copysign(0.5 − 0.25ε, x))` | A | 659 | **all 0** | **all 0** (x86_64 PE too) |
| `ldexp` | **`_tm_ldexp`** (M2): Rust libm's scalbn | A | 306 | **all 0** | **all 0** (x86_64 PE too) |
| `sin` | `f64_sin` | C | 710 | all 0 | all 0 |
| `cos` | `f64_cos` | C | 710 | all 0 | all 0 |
| `tan` | `ganita_f64_tan` | C | 710 | all 0 | all 0 |
| `exp` | `f64_exp` | C | 710 | 650 × 0, **60 × 1** | all 0 |
| `log` | `f64_ln` | C | 874 | 864 × 0, 10 × 1 | 870 × 0, 4 × 1 |
| `sinh` | `ganita_f64_sinh` | C | 510 | 495 × 0, 15 × 1 | 509 × 0, 1 × 1 |
| `pow` | `ganita_f64_pow` | C | 1,220 | 1,206 × 0, 13 × 1, 1 accepted | same |
| `hypot` | `ganita_f64_hypot` (measured in M0) | C | 600 | 503 × 0, 96 × 1, 1 × 2 | same |
| `hypot` | **`_tm_hypot`** (from M1) | **A** | 600 | **all 0** | **all 0** (x86_64 PE too) |
| `sin` | **`_tm_sin`** (from M4): Rust libm's sin, k_sin, k_cos and rem_pio2 transcribed | **A** | 710 | **all 0** | **all 0** (x86_64 PE too) |
| `cos` | **`_tm_cos`** (from M4) | **A** | 710 | **all 0** | **all 0** (x86_64 PE too) |
| `log` | **`_tm_log`** (from M4): Rust libm's log (musl's form) transcribed | **A** | 874 | **all 0** | **all 0** (x86_64 PE too) |

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

**Done in M1:** `src/numerics.cyr` has `_tm_hypot`. The harness's `libm::hypot` rows go to it and are class A,
600/600 bit-identical on x86_64 ELF, x86_64 PE and aarch64. With it, every four-momentum golden row is 0 ulp on all
three targets, including the 638 ill-conditioned `invariant_mass` rows.

## Class D: Bateman (M2)

`bateman_chain` is a scaling-and-squaring matrix exponential: exp and sinh on the diagonals, then up to about 70
squarings for the U-238 series at 4.463e9 y. Over its 28 golden rows (183 populations), on 2026-11-08:

| Target | ulp histogram | max |
|---|---|---|
| x86_64 ELF and PE | 125 × 0, 29 × 1, 18 × 2, 11 × 3–16 | 9 ulp (about 2e-15 relative) |
| aarch64 | all 0 | 0 |

The difference is `f64_exp`: x87 on x86_64, and the fdlibm polyfill on aarch64, which is the same algorithm Rust
libm uses. Class D's relative 1e-12 holds with about three decades to spare. `_tm_ldexp`, the squaring scale, is
Rust libm's scalbn transcribed, and is bit-identical (306/306 `libm::ldexp` rows, class A, on every target).

## Atomic (M3)

Over the 43,309 `atomic.tsv` rows on 2026-11-08, every class A and B row is 0 ulp on all three targets: the
configurations and formats, the Rydberg and reduced-mass lines, fine structure, the series, Zeeman and Stark,
vacuum polarization, Dirac, hyperfine, the g-factors and Breit. The class-C and class-D functions:

| Function | Class (K) | Rows | x86_64 ELF and PE | aarch64 |
|---|---|---|---|---|
| `lamb_shift_nlj_ev` | C (10) | 2,596 | 1,906 × 0, 12 × 1, 2 × 2 | 1,918 × 0, 2 × 1 |
| `lamb_shift_ev` | C (10) | 2,520 | 2,506 × 0, 12 × 1, 2 × 2 | 2,518 × 0, 2 × 1 |
| `radial_wavefunction` | C (6) | 9,249 | 8,298 × 0, 768 × 1, 129 × 2 | all 0 |
| `radial_probability_density` | C (18) | 9,249 | 8,319 × 0, 189 × 1, 356 × 2, 331 × 3–7 | all 0 |
| `einstein_a_coefficient` | D | 1,443 | all 0 | all 0 |
| `einstein_b_coefficient` | D | 1,443 | all 0 | all 0 |

(The histograms count float outputs; Err rows have none.)

- **The targets differ only through exp and ln**, as for the primitives: aarch64's fdlibm exp makes the radial
  functions bit-identical there, and its ln polyfill leaves 2 one-ulp Lamb rows against 14 on x86 (x87).
- **The radial dipole integral is bit-identical on this grid.** Einstein A and B are class D: a double sum of pow
  terms whose alternating signs cancel. Every row here matches exactly, out to n_u = 60, with `ganita_f64_pow` for the
  (2/n)^i factors. The M6 optics rows found the exception: of the 4,215 transitions `spectral_series_lines` uses
  (Z ∈ {1, 2, 26}, n_l ≤ 5, n_u ≤ 60), three differ, all 49p → 1s, by 55–96 ulp (about 2 × 10⁻¹⁴ relative). The
  49-term sum cancels, and pow is the only primitive on that path that isn't exact. Class D's 1e-12 holds by two
  orders of magnitude.
- **Dirac is class B, not C.** `dirac_energy_mev` and its two callers use only sqrt and arithmetic, so their
  2,856 rows are B and bit-identical. That tightens the plan's class.

## Reactions and scattering (M4): exact sin, cos and log

The port-parity rule measured S by perturbation (each sin, cos, tan and ln result nudged ±1 ulp, one call at a time,
over every golden row) before the first comparison. Twelve of the sixteen class-C candidates came out past K = 64:

| Function | S (worst row) | Where |
|---|---|---|
| `average_lethargy_gain` | 4.3 × 10⁹ | ξ = 1 + ((A−1)²/2A) ln((A−1)/(A+1)) is 1 + (≈ −1) at large A |
| `collisions_to_thermalize` | 211 | through ξ |
| `mott_correction_factor` | 4.4 × 10¹⁸ | 1 − β² sin²(θ/2) at β = 1, θ = π |
| `mott_electron_differential`, `mott_electron_with_form_factor` | 1.5 × 10⁹ | the same factor at 10 GeV, θ = π |
| `nuclear_form_factor_uniform` | 17,337 | near its first zero, qR = 4.4934 |
| `compton_energy_ratio` | 4.6 × 10¹⁸ | the ∞ · 0 edge (E_γ = ∞, θ = 0); finite rows near θ = 0 are past 64 too, through 1 − cos θ |
| `klein_nishina_differential` | 77,700 | 10 GeV at θ = 0, through 1 − cos θ |
| `klein_nishina_total` | 12,303 | the closed form just above γ = 1e-2 |
| `partial_wave_differential` | 187 | eleven phase shifts at θ = π |
| `pair_production_cross_section` | 73 | just above threshold, Z = 82 |
| `partial_wave_cross_section` | 3, but 11 sin calls: K = 68 | |

Step 3 of the rule applies, as for hypot: **the port reproduces the primitives bit for bit.** `_tm_sin`, `_tm_cos` and
`_tm_log` (src/numerics.cyr) transcribe libm 0.2.16's sin.rs, cos.rs, k_sin.rs, k_cos.rs, rem_pio2.rs and log.rs
operation for operation. Two differences from the stdlib made that necessary:
- **`f64_ln`** is x87 `fyl2x` on x86_64, and on aarch64 a port of fdlibm's branchy e_log.c. Rust's log is musl's
  single-expression form. A sweep found 16,248 of 680,010 inputs differing on x86_64 and 2,712 on aarch64.
- **`f64_sin` and `f64_cos`** share Rust's kernels and constants, but their argument reduction always takes the
  "medium" path, where Rust has direct reductions for |x| ≤ 9π/4. The two differ only when |x − kπ/2| is in about
  [2⁻²⁰, 2⁻¹⁶)·|x| and the result sits within a hair of a rounding tie. The sweep found none, so this is closed by
  construction rather than by sample.

Arguments at or past 2²⁰·π/2 go to the stdlib's Payne–Hanek reduction (`_f64_rem_pio2`, FreeBSD k_rem_pio2.c, as Rust's
rem_pio2_large is). The sweep covers 40,000 of them, and libm.tsv has 29.

**The sweep.** 680,010 inputs aimed at the risky regions: random bit patterns, (π/4, 9π/4), the doubles near kπ/2 for
k ≤ 12 at relative offsets 2⁻¹⁰ … 2⁻⁵², medium and huge magnitudes, log near 1 and across the subnormals. `_tm_sin`,
`_tm_cos` and `_tm_log` match Rust libm on **every input on x86_64 ELF, x86_64 PE and aarch64**. The cost: `_tm_sin`
takes 48 ns against `f64_sin`'s 37 ns, and `_tm_log` 54 ns against x87's 38 ns.

**Results.** With the exact primitives, every reaction and scattering function except one is class B, and all 18,669
rows are **0 ulp on all three targets**. `rutherford_total_above_angle` keeps ganita's tan, with S = 4 and K = 10 (below),
and is also 0 ulp. The `libm::sin`, `cos` and `log` rows move to the new primitives at class A, a tightening, and are
exact on every target (table above). The M0 rows for `f64_sin`, `f64_cos` and `f64_ln` keep their record of what the
builtins do.

## Optics (M6)

The 1,039 `optics.tsv` rows come from the Rust oracle linked with the Rust prakash 1.2.0; the port runs against the
Cyrius prakash 2.15.2, replayed by tests/optics_parity.tcyr. All pass on all three targets, on the first comparison
(made against 2.15.1 in 2.0.0). 2.0.1 moved to 2.15.2, whose 21 functions that `src/optics.cyr` reaches are
byte-identical to 2.15.1's, so the table holds:

| Function | Class (K) | Rows | x86_64 ELF and PE | aarch64 |
|---|---|---|---|---|
| `hydrogen_balmer_series` | B | 1 (12 values) | all 0 | all 0 |
| `spectral_series_lines` | D | 168 | 2,247 × 0, 3 × 55–95 | the same |
| `lines_to_spd` | C (38) | 32 (4,391 values) | 4,268 × 0, 113 × 1, 10 × 2 | **all 0** |
| `line_to_rgb` | B | 838 | all 0 | all 0 |

- **prakash agrees with its Rust self.** `line_to_rgb` is prakash's piecewise CIE fit, so 838 bit-identical rows
  across the 1.2.0 → 2.15.2 port boundary, every branch edge included, say the two prakash lines compute it the same
  way.
- **The SPDs differ only through exp**: x87 `f64_exp` on x86_64, fdlibm on aarch64 as in Rust. A sum of positive
  Gaussians is well-conditioned, so 2 ulp is the most it moves.
- **The three D outliers are Einstein A at 49p → 1s** (above, M3), through the series' A-value intensities.

## Chrono bridge (2.0.1)

The 6,588 `bridge.tsv` rows are the first bridge's (port plan §6). All are **0 ulp on all three targets**, on the
first comparison:

| Function | Class | Rows | x86_64 ELF and PE | aarch64 |
|---|---|---|---|---|
| `utc_date_to_tai_offset` | A | 1,170 | all equal | all equal |
| `tai_to_utc_seconds` | B | 2,709 | all 0 | all 0 |
| `utc_to_tai_seconds` | B | 2,709 | all 0 | all 0 |

- The two conversions are one IEEE add or subtract of ΔAT, so class B and 0 ulp were the expectation. The grid crosses
  the float specials and the 2⁵³ neighbourhood, where adding ΔAT rounds, with every month where ΔAT changes and the
  month before it.
- `validate.py` recomputes every conversion row in Python and finds the same bits. It checks the offsets against
  `timekeeping.tsv`'s `leap_seconds_at` rows, and against the IERS shape: 0 before 1972, 10 at 1972-01, +1 s steps,
  37 from 2017-01.

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
| `Nucleus::atomic_mass_amu` | M2 | **0**: the two `pow` results (Z^2.39 and Z^5.35) nudged ±1 ulp on all 2,053 golden rows, with `nuclear_mass` taken from the golden corpus. The electron binding is about 10⁻⁶ of the mass, so a 1-ulp change never moves the output | 2 | **2** |
| `sagnac_correction_ns` | M5 | **2**: sin(φ) nudged ±1 ulp on all 21 golden rows, through four products and a quotient | 1 | **6** |
| `remaining_fraction` | M2 | **1**: the output is the pow result itself, `pow(0.5, t / T½)` (t / T½ is a correctly rounded division on both sides) | 1 | **4** |
| `lamb_shift_nlj_ev` | M3 | **4**: ln(1/(Zα)²) nudged ±1 ulp on every s-state golden row (Z = 1–118 and u32::MAX, n ≤ 4). The l > 0 rows have no ln (their coefficient is the tabulated Bethe logarithm), and Z = 0 is 0 · ln(∞) = NaN whatever ln returns | 1 | **10** |
| `lamb_shift_ev` | M3 | **4**: an s state is `lamb_shift_nlj_ev(z, n, 0, 1)` itself; l > 0 averages two ln-free values | 1 | **10** |
| `radial_wavefunction` | M3 | **2**: exp(−ρ/2) nudged ±1 ulp on every golden row; R = (N·P)·e, so a 1-ulp change of e moves R by at most 2 ulp across a binade. 26 edge rows whose exp result is 0 or ∞ (r ≥ 1e300, ±∞, NaN) are exact on both sides and excluded | 1 | **6** |
| `radial_probability_density` | M3 | **8**: as `radial_wavefunction`, through R² r² | 1 | **18** |
| `rutherford_total_above_angle` | M4 | **4**: ganita's tan(θ_min/2) nudged ±1 ulp on all 990 golden rows; the output goes as 1/tan² | 1 | **10** |
| `lines_to_spd` | M6 | **3**: each exp result nudged ±1 ulp on every sample of all 32 golden SPDs (a 1-line, 2.5-intensity SPD is the worst); an exp that underflows to 0 is exact on both sides | 6 (one exp per line, at most 6 lines) | **38** |

The Myers–Swiatecki shell term turned out to be class B, not C. It uses only cbrt (correctly rounded on both
sides) and arithmetic, and its 14,160 rows are bit-identical.

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
