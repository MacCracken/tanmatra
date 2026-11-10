# Reference validation

The golden corpus checks the port against Rust 1.3.0 ([`port-parity.md`](port-parity.md)). This file covers the other
check: tanmatra against evaluated primary data, where no Rust row exists or where Rust was never checked. Each section
gives the reference, the functions compared, and each bound with its derivation.

## The rule for a bound

1. Name what the function leaves out, and take the leading omitted term analytically, as its size for the inputs
   compared.
2. The bound is **twice that size**. List every further omitted term with its size, and check that together they fit
   well inside the second half.
3. Write the bound and its derivation here **before the first comparison**. A residual past the bound is a finding.
   Either the code has a defect, which is fixed (and needs an ADR when the fix changes what Rust 1.3.0 returned), or
   the derivation missed a term. A missed term is added here with its physics, and the bound is never fitted to a
   residual.
4. After the run, record the worst residual as a fraction of its bound.
5. Record what each bound catches, so a bound loose enough to pass anything is visible.

The reference's stated uncertainty is one more term in the second half's budget. Each test asserts that it is at most
**a tenth of the bound** before it compares. That leaves room for the next omitted terms, and keeps a residual a
measurement of tanmatra, not of the reference.

## He⁺ and Li²⁺ (2.0.2)

**Reference.** NIST ASD v5.12 energy levels of He II (⁴He) and Li III (⁷Li). That is 148 levels per ion: every (n, l,
j) for n ≤ 11, and for n = 12–20 the nS½ and nP½ levels and the one with the highest l. They come from Yerokhin &
Shabaev 2015 (n ≤ 2 and the limit) and Erickson 1977, scaled by NIST (n ≥ 3). The snapshot and its provenance are in
[`tests/reference/nist_asd/`](../../tests/reference/nist_asd/README.md). `scripts/reference/nist_levels.py` generates
`tests/nist_levels_data.cyr` from it. The test is `tests/hydrogenic_nist.tcyr`, at (Z, A) = (2, 4) and (3, 7).

**Notation.** R = R∞ (`TM_RYDBERG`, m⁻¹), α = `TM_FINE_STRUCTURE`, and μ = `reduced_mass_factor(z, a)`. Wavenumbers
are in m⁻¹, and NIST's cm⁻¹ are multiplied by 100. A level (n, l, j) has NIST energy E above the ground level 1s ²S½.

### 1. Bohr lines: `spectral_line_vacuum_nm` and the named series

**Compared.** Every pair of NIST levels (n₁, l₁, j₁) and (n₂, l₂, j₂) with n₁ < n₂ ≤ 20, which is 10,026 pairs per ion.
The test compares σ = 10⁹ / `spectral_line_vacuum_nm(z, a, n₁, n₂)` against E₂ − E₁. The named series
(`lyman_series` … `pfund_series`, infinitely heavy nucleus) divided by μ must equal `spectral_line_vacuum_nm` bit for
bit for every n ≤ 20, so the same comparison holds them.

**Leading omitted term: the fine structure.** To first order, level (n, j) is more bound than Bohr's by

  F(n, j) = R Z² (αZ)² / n⁴ × (n / (j + ½) − ¾),

which is positive for every j ≤ n − ½. Each component of the line moves by F(n₁, j₁) − F(n₂, j₂), so its size is
below F(n₁, j₁) + F(n₂, j₂).

**Bound.** B₁ = 2 [F(n₁, j₁) + F(n₂, j₂)], per pair, from each level's own j. μ is left out of F, which makes the bound
larger by at most 1/μ − 1 ≤ 1.4 × 10⁻⁴ of itself.

**Next terms**, each as a fraction of the F of the level it shifts:
- the one-loop QED (Lamb) shift: (α/π)(Zα)⁴ m_ec² / n³ × [4/3 ln(Zα)⁻² − 4/3 ln k₀ + 10/9 − 4/15 + Zα A₅₀] for s
  levels. It is 0.154 at 1S and 0.063 at 2S for He⁺ (0.135 and 0.056 for Li²⁺), and falls with n. For l > 0 it is
  at most about 2(α/π) ≈ 0.005;
- the higher-order Dirac terms, ~(Zα)² ≤ 4.8 × 10⁻⁴;
- recoil beyond the reduced mass, ~m_e/M ≤ 1.4 × 10⁻⁴;
- the finite nuclear size, s levels only: (16/3)(r_N/ƛ_C)² at 1S, 1.0 × 10⁻⁴ for ⁴He and 2.1 × 10⁻⁴ for ⁷Li.

Together they are at most 0.16 of F, so a correct function stays below 0.58 B₁. NIST's stated uncertainties, summed
over the two levels, are at most 6.0 × 10⁻⁵ of B₁.

**What it catches.** Relative to the line, B₁ is as small as 1.3 × 10⁻⁶ (He⁺) and 2.9 × 10⁻⁶ (Li²⁺), between the
high-l levels at large n. Any reduced-mass error at least that size fails:
- no reduced mass (μ = 1): 1.4 × 10⁻⁴ and 7.8 × 10⁻⁵;
- the proton's mass in place of the nucleus's: 4.1 × 10⁻⁴ and 4.7 × 10⁻⁴;
- ³He or ⁶Li in place of ⁴He or ⁷Li: 4.5 × 10⁻⁵ and 1.3 × 10⁻⁵.

A wrong power of Z is off by a factor of order 1.

### 2. Fine-structure splittings: `spectral_line_fine_nm`

**Compared.** For n = 2–11 and l = 1 … n − 1 (55 splittings per ion), the test compares
σ = 10⁹ / `spectral_line_fine_nm(z, n, 2l − 1, n, 2l + 1)` against E(n, l, l + ½) − E(n, l, l − ½).

**The formula.** It is first-order and infinite-mass: ΔF = R Z² (αZ)² / (n³ l (l + 1)).

**Leading omitted term: the electron's anomalous magnetic moment.** This is the κ term of the one-loop QED shift,
−(α/π)(Zα)⁴ m_ec² / n³ × 1/(2κ(2l + 1)). Between j = l ± ½ it adds exactly (α/π) ΔF, for every n and l. The
j-independent part of the Lamb shift (the Bethe logarithm) cancels in the splitting.

**Bound.** B₂ = 2 (α/π) ΔF, with ΔF from R and α, not from the function under test.

**Next terms**, as fractions of ΔF:
- the higher-order Dirac terms, ~(Zα)² ≤ 4.8 × 10⁻⁴;
- the reduced mass, ≤ 1.4 × 10⁻⁴;
- the two-loop terms, ~(α/π)² ≈ 5 × 10⁻⁶.

Together they are under 0.3 (α/π), so a correct function stays below 0.65 B₂. NIST's stated uncertainties, summed
over the two levels, are at most 0.058 of B₂. The sum overstates the splitting's uncertainty: NIST gives each level's
as a separation from the limit, and the shared part cancels in a splitting.

**What it catches.** B₂ is 4.6 × 10⁻³ of the splitting. A wrong power of Z, n or l, or a wrong j dependence, is off
by order 1.

### 3. The ⁴He⁺ reduced mass: `reduced_mass_factor(2, 4)`

**Compared.** The test compares μ against 1 / (1 + m_e/m_α), with CODATA 2022 m_α/m_e = 7294.299 541 71(17)
(https://physics.nist.gov/cgi-bin/cuu/Value?malsme).

**Leading omitted term.** The function takes M as the AME2020 ⁴He atomic mass less two electron masses. That drops
the electrons' binding energy, 79.005 eV (He I 24.587 eV plus He II 54.418 eV, NIST ASD), which is 2.12 × 10⁻⁸ of
m_α c². μ moves by (m_e/M) × 2.12 × 10⁻⁸ = 2.9 × 10⁻¹² of itself.

**Bound.** B₃ = 5.8 × 10⁻¹² relative.

**Next terms.** The CODATA uncertainty moves μ by 3 × 10⁻¹⁵ and AME2020's ⁴He uncertainty by 5 × 10⁻¹⁵. Rounding is
~10⁻¹⁶.

**What it catches.** Any wrong mass: ³He is 4.5 × 10⁻⁵ off, and a dropped Z m_e term is 3.8 × 10⁻⁸ off.

There is no CODATA ⁷Li mass. For Li²⁺, §1's high-n pairs hold μ to better than 3 × 10⁻⁶.

### Results (2026-11-09)

Every comparison passes on x86_64 ELF, x86_64 PE and aarch64 (61,948 assertions, the same ratios on all three). The
worst residual of each bound is the term the derivation named:

| Bound | He⁺ | Li²⁺ | The worst residual is |
|---|---|---|---|
| B₁ (10,026 pairs per ion) | 0.5004 | 0.5005 | 2P½ to the n = 20 circular level: 2P½'s fine structure plus its own one-loop shift, 0.4 (α/π) of it. The circular level's F is 2 × 10⁻⁵ of 2P½'s, so the ratio sits at ½ |
| B₂ (55 splittings per ion) | 0.4980 | 0.5466 | the anomalous moment. Every splitting's residual is 0.93–1.00 (He⁺) and 0.96–1.09 (Li²⁺) times (α/π) ΔF, and the spread is NIST's uncertainty and the (Zα)² terms |
| B₃ | 0.5008 | — | the dropped 79.005 eV electron binding, 2.9 × 10⁻¹², as derived |

The first run set the uncertainty precondition at a hundredth of the bound, a number not taken from the budget. 32 of
the 55 splittings per ion failed it, though none failed B₂. The rule now sets the precondition at a tenth, as part of
the second half's budget. B₁, B₂ and B₃ are unchanged from their derivations.
