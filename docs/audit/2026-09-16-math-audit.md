# Audit — 2026-09-16 (v1.2.1): math and data correctness before the Cyrius port

**Scope.** Every formula and every tabulated value in `src/`. Findings come from
running code, not from reading it: a scratch harness calls the crate as a path
dependency. Its output was compared against:

- **Independent computation:** exact normalization integrals, analytic limits,
  and 60–120-digit `Decimal` evaluations.
- **Machine-readable primary data:**
  - CODATA 2022 complete listing
  - AME2020 `mass_1.mas20`
  - NUBASE2020 `nubase_4.mas20`, IAEA LiveChart, DDEP
  - ENDF/B-VIII.0 evaluation files, with 0 K cross sections rebuilt from their
    resonance parameters
  - Sears 1992, IAEA Standards 2017
  - Angeli & Marinova 2013; Stone INDC(NDS)-0658/0794/0833
  - Hardy & Towner 2009/2015/2020 survey texts
  - NIST ASD ionization energies, Andersen–Haugen–Hotop 1999 plus newer measurements
  - IERS `Leap_Second.dat`, IANA `leap-seconds.list`, USNO `tai-utc.dat`
  - BIPM standard-frequency API, IERS Conventions 2010
  - PDG 2024 MC table

Four data-table sweeps (atomic, decay, reaction/nuclear, timekeeping/particle) ran as
independent agents. **Every high-severity agent finding was re-derived from the
primary file before it was allowed into this report.** No finding was refuted on
re-derivation.

**Severity.**
- **high**: a physically wrong result, or the wrong physical quantity, in normal use.
- **medium**: wrong beyond the reference uncertainty in a way that matters, or
  wrong only in part of the input domain, or documentation that leads to wrong use.
- **low**: numerically negligible, provenance only, validation, or robustness.

## Project state (catch-up)

- **Version:**
  - `Cargo.toml`/`VERSION` say **1.2.1**, but `CHANGELOG.md` stops at **1.2.0**.
  - The 1.2.0 release shipped roadmap v1.2–v1.5 (nuclear data tables, advanced
    scattering, relativistic quantum, frequency standards and time).
  - Open: v1.1 (wavefunctions for n > 3, Y_lm, radial integrals, oscillator
    strengths) and v1.6 (simulation support).
- **Health:**
  - `fmt`, `clippy -D warnings` and `doc -D warnings` are clean.
  - **451 tests pass** (429 unit, 20 integration, 2 doc) with all features; 439 with
    no default features.
  - `cargo audit` **fails** on a dev-only dependency: crossbeam-epoch 0.9.18 via
    criterion 0.5 (RUSTSEC-2026-0204). `cargo deny` fails on the same advisory.
- **Benchmarks:** 11, all from 1.0.0. Nothing added in 1.1 or 1.2 is benchmarked
  (timekeeping, the new scattering, Dirac/QED, data lookups).
- **Stale documentation:**
  - `docs/guides/testing.md` says 191 unit and 212 total tests, and 91.7% coverage.
  - `docs/architecture/overview.md` does not list `timekeeping`, `bridge` or
    `integration`.
  - `roadmap.md` leaves the shipped bridges and soorat sections unchecked, and
    checks TCG/TCB conversions that do not exist (D6).
  - The feature table in `lib.rs` omits `optics` and `soorat-compat`.

## Headline

**Formulas: 36 confirmed** — 14 high, 15 medium, 7 low.
**Data tables: 62 confirmed rows** — 22 high, 19 medium, 21 low. Some rows group
several elements or nuclides.

1. **Scattering functions:**
   - `distance_of_closest_approach` is 2× too small.
   - `rutherford_total_above_angle` is 4× too small.
   - `pair_production_cross_section` goes **negative** from threshold to 3.4 MeV.
2. **Atomic formulas:**
   - One hydrogen wavefunction (R₃₁) integrates to **36** instead of 1.
   - The Einstein A coefficient is off by factors of 0.025 to 21 and misses its own
     calibration point by 3×.
   - Vacuum polarization is 25% high.
   - The Breit term exceeds the entire measured two-electron energy of He-like
     uranium 15-fold.
3. **The "Strutinsky shell correction"** makes binding energies worse for **98.6%**
   of 2484 measured nuclides. ADR-005 claims the opposite.
4. **Shell model:** every odd-odd nucleus gets an impossible half-integer spin
   (352 of 352).
5. **`bateman_chain`:**
   - It silently drops terms: for equal decay constants the answer is wrong and does
     not conserve atoms.
   - For long-lived chains the result depends on whether time is in seconds or years.
   - At short times, trace daughters come out 100× wrong.
6. **Provenance:** almost no table is the edition it cites.
   - "CODATA 2022" is 2018.
   - Only 52 of 114 half-lives are NUBASE2020.
   - The thermal cross sections are Sears 1992 **bound-atom** values; the H-1
     scattering entry is 4× the free-atom value.
   - The superallowed data is the 2009/2015 surveys, not 2020.
   - The electron affinity table has ~30 wrong or superseded entries, including
     measured bound anions reported as unbound.
7. **The test suite cannot catch this class of error.** Tolerances are wide (Einstein
   A accepts 1e7–1e10). Some tests re-derive the formula under test (Breit). Others
   pin the wrong values:
   - closest approach 22.8 fm
   - Fe-56 mass excess −60601 keV
   - C-14 at 5730 ± 1 y
   - U-235 "capture" resonance integral 275 b
   - Sr frequency ± 10 Hz

   None of the 451 tests fails on any finding below.

---

## Part 1 — Formulas

### High

| # | Where | Defect | Reproduction |
|---|---|---|---|
| F1 | `scattering.rs:112` `distance_of_closest_approach` | Returns Z₁Z₂e²/(2T); the head-on turning point is Z₁Z₂e²/T. Test `closest_approach_gold_alpha` pins the wrong 22.8 fm | 5 MeV α on Au: code 22.75 fm, correct 45.50 fm |
| F2 | `scattering.rs:83` `rutherford_total_above_angle` | Returns π a²/tan²(θ/2) with a = Z₁Z₂e²/4E. The correct value is 4π a²/tan²(θ/2) = π (Z₁Z₂e²/2E)² cot²(θ/2) | numeric ∫(dσ/dΩ)dΩ over θ > 0.1, divided by code: **4.0000** |
| F3 | `scattering.rs:461` `pair_production_cross_section` | The unscreened high-energy Bethe–Heitler asymptote is used from threshold upward. It is **negative** for 1.022 < E < 3.42 MeV, has no screening at high E (~50% high for Pb at 1 GeV), and has no electron-field term | Pb at 1.05 MeV: −14.3 b; at 2 MeV: −6.5 b |
| F4 | `atomic.rs:829` `radial_wavefunction(n=3, l=1)` | Uses (6 − Zr) where the normalized form has (1 − Zr/6): R₃₁ is 6× too large | ∫R²r²dr = **36.000000** (all other n, l give 1.000000, for Z = 1 and 3) |
| F5 | `atomic.rs:937` `einstein_a_coefficient` | Heuristic A_Lyα·Z⁴·(ν/ν_Lyα)³·l_max/(2l_u+1), not the dipole formula. The angular factor breaks its own calibration point, and scaling goes as Z¹⁰ because ν already carries Z² | H Lyα 2.09e8 vs NIST 6.2649e8 (×0.33); Lyβ ×2.08; 3p→2s ×0.059; 3s→2p ×0.63; 3d→2p ×0.025; He⁺ Lyα ×21.4 |
| F6 | `atomic.rs:754` `vacuum_polarization_ev` | Coefficient α/(3π). The Uehling s-state shift is −(4/15)(α/π)(Zα)⁴mc²/n³, so the code is 25% high | H 2S: code −33.91 MHz; Uehling −27.13 MHz |
| F7 | `atomic.rs:1542` `breit_interaction_ev` | −(2/3)α²Z⁴·Ry has the Z-scaling of a one-body relativistic correction; a two-electron interaction scales as α²Z³ at leading order. Unsourced; its test re-derives the same formula | U⁹⁰⁺: code **−34.6 keV**, while the *entire* measured two-electron contribution is 2248 ± 9 eV (Gumberidze et al., PRL 92, 203004, 2004) |
| F8 | `nucleus.rs:277` `binding_energy_shell_corrected` | A sum of positive Gaussians with unsourced peak energies, labelled "Strutinsky" (it is not). It can only add binding, never remove it, even though the docs say mid-shell corrections are negative | vs AME2020, 2484 nuclides with A ≥ 16: RMS **10.24 → 14.58 MeV**; improves 36 of 2484; doubly-magic mean error **+0.69 → +7.79 MeV** |
| F9 | `nucleus.rs:673` `ground_state_spin_parity` | Odd-odd nuclei return the odd proton's j: a **half-integer spin for even A**, and a parity that ignores the neutron | **352 / 352** firm NUBASE2020 odd-odd J^π wrong, all half-integer (N-14 gives 1/2⁻; measured 1⁺) |
| F10 | `decay.rs:404` `bateman_chain` | Terms are silently skipped for degenerate λ. The absolute thresholds `1e-30`/`1e-100` make the result depend on the time unit. The alternating sum cancels catastrophically at short times | λ = [0.1, 0.1, 0], t = 10: code **[0.368, 0.0, 1.0]**, exact [0.368, 0.368, 0.264] (atoms not conserved). λ = 1…9 × 1e-13 s⁻¹, t = 1e13 s: last member 0.368 in seconds, 0.009378 in years (exact 0.009378). U-238 series in SI at t = 1 y: Po-210 131× high, Pb-210 13× high, Pb-206 = 0 (all correct at t ≥ 10³ y) |
| F11 | `decay.rs:159` `decay_chain` | `Nucleus` cannot tell an isomer from its ground state. IT maps a nucleus to itself and the lookup finds the isomer again, so the chain loops | Ta-180, Hf-178, Am-242: `IsomericTransition` repeated on the same nucleus until `max_steps` |
| F12 | `bridge.rs:42` `atomic_number_to_valence` | `_ => 2` for every Z ≥ 19 (and Z = 0) | K → 2 (should be 1), Ga → 2 (3), Br → 2 (7), Kr → 2 (8), I → 2 (7), Cs → 2 (1) |
| F13 | `timekeeping.rs:379` `AtomicInstant::add_seconds` | Round-trips through f64 seconds since 1958 (ulp ≈ 4.8e-7 s today). This destroys the nanosecond field the type exists for, and is not even the identity for dt = 0 | at 2.1e9 s, both `add_seconds(0.0)` and `add_seconds(1e-9)` move the instant **−73 ns** |
| F14 | `timekeeping.rs:479` `sagnac_correction_ns` | 4ΩA·sinφ/c² is the ring-laser formula (closed loop, two beams, horizontal). The doc describes a one-way signal, which is 2ΩA_E/c², and sinφ re-projects an area already described as projected | equatorial loop: code 0.000 ns; one-way Sagnac delay 207.39 ns (Ashby 2006) |

### Medium

| # | Where | Defect | Reproduction |
|---|---|---|---|
| F15 | `atomic.rs:407` `spectral_line_nm` (and every series, `optics.rs`, soorat) | Uses R∞ with no reduced-mass factor 1/(1+mₑ/M): 5.4e-4 relative for H. The threat model claims "< 0.01%, exact for hydrogen-like". The docs and tests quote 656.3 nm, which is the **air** wavelength | Hα: code 656.112 nm; reduced-mass vacuum 656.470 nm |
| F16 | `atomic.rs:461, 962, 1015` | Hard-coded `13.6` eV (4.2e-4 low) in `hydrogen_level_energy_ev` and both Einstein coefficients, although `RYDBERG_EV` exists in the same file. The truncation is ~8× the fine-structure effect the function models | E(1s½) −13.60018 vs −13.60587 eV; `spectral_line_fine_nm` moves Hα by 0.27 nm (actual fine structure ≈ 0.016 nm) |
| F17 | `atomic.rs:1446` `hyperfine_splitting_ev` | Missing the (I + ½) factor, so only I = ½ is right | D 1s: code 218.16 MHz, measured 327.384 MHz; ×1.5 gives 327.24 |
| F18 | `atomic.rs:716` `lamb_shift_ev` | Pure Z⁴/n³ scaling of the H 2S value (no ln(Zα)⁻², no F(Zα)). `0.1/l` for l > 0 is invented and has the wrong sign. `4.3725e-6` ≠ h·1057.845 MHz (4.37488e-6 eV) | H 1S +3.5%; He⁺ 2S +20%; U⁹¹⁺ 1s 2506 eV vs measured 460 eV; H 2P +105.7 MHz vs ≈ −12.8 MHz |
| F19 | `nucleus.rs:316` `atomic_mass_amu` | Returns the **nuclear** mass/u; the atomic mass adds Z·mₑ | Fe-56: 55.91359 u (the atomic mass would add 0.01426 u) |
| F20 | `nucleus.rs:176` SEMF coefficients | A textbook set fitted with a Z² Coulomb term, used here with Z(Z−1). Systematic +9.2 MeV overbinding; ADR-005/threat model claim ~1–2% | The same 5-term form least-squares fitted to AME2020 gives RMS **3.31 MeV** vs 10.24 (a_v 15.414, a_s 16.860, a_c 0.6952, a_a 22.497, a_p 12.03) |
| F21 | `nucleus.rs:585` shell level order | The 82–126 block puts 1i13/2 last; the 50–82 block puts 1g7/2 first | odd-A J^π agrees with **38.4%** (326/849) of firm NUBASE2020 assignments; Pb-207 gives 13/2⁺ (1/2⁻), I-127 7/2⁺ (5/2⁺), Sn-117 3/2⁺ (1/2⁺) |
| F22 | `nucleus.rs:776` `corrected_ft_value(_ft)` | Ignores its argument and returns a constant — a stub documented as applying δ_R′, δ_NS, δ_C | `corrected_ft_value(1.0)` = 3072.27 |
| F23 | `scattering.rs:305` `mott_electron_differential` | Ultra-relativistic limit (E ≈ pc, β = 1), documented as taking kinetic energy "(relativistic)" | 1 MeV kinetic e⁻ on Au at θ = 1: 1.73× the point-nucleus Mott value |
| F24 | `scattering.rs:69` `mott_correction_factor` | (1 − β² sin²θ/2) is the relativistic spin-½ Mott factor; the docs describe identical-particle exchange, which is a different formula | docs vs physics |
| F25 | `integration/soorat.rs:35` `hydrogen_slice` | No Laguerre polynomial (no radial nodes), no Y_lm (no shape), ρ scaled differently from `radial_wavefunction` | 2s slice has no node at r = 2a₀; 2p density is identical on the x and z axes |
| F26 | `bridge.rs:111` `nuclear_spin_to_hyperfine_scale` | Returns the multiplicity 2I+1 as a "splitting factor"; the splitting goes as μ_I(I+½)/I | unphysical |
| F27 | `atomic.rs:1516` `anomalous_zeeman_splitting_ev` | Implements (m_l + g_e m_s)μ_B B, the Paschen–Back (strong-field) formula; the doc header says m_j·g_e·μ_B·B | naming/docs |
| F28 | `scattering.rs:246` `born_screened_coulomb` | Masses guessed as A ≈ 2Z (proton mass for Z = 1), so electron projectiles are impossible. The doc formula has a factor-4 error and a stray ħ; the code itself is right | Au: A = 158 used vs 197 |
| F29 | `timekeeping.rs:252` `leap_seconds_at` (and the three bridge UTC functions) | Returns 0 before 1972 | TAI−UTC was 1.42 s on 1961-01-01 and 9.892 s at the end of 1971 |

### Low

| # | Where | Defect |
|---|---|---|
| F30 | `scattering.rs:423` `klein_nishina_total` | Cancellation for 1e-6 ≤ γ ≲ 1e-4: 3.8% error at γ = 1e-5. Needs a series branch |
| F31 | `scattering.rs:331` `nuclear_form_factor_uniform` | Cancellation near its 1e-6 cutoff: 7.8e-5 error at qR = 1e-6. Needs a series branch below ~1e-2 |
| F32 | `atomic.rs` validators | `hydrogen_level_energy_ev` accepts integer j; `dirac_energy_mev` accepts j > n−½; `lande_g_factor` and `stark_shift_hydrogen_ev` accept impossible (l, j) and (n, k); `spectral_line_nm(z=0)` returns `Ok(inf)` |
| F33 | `relativity.rs` | `lorentz_gamma(−1.5)` = NaN (checks β ≥ 1, not \|β\| ≥ 1); `velocity_addition(1, −1)` = NaN; `invariant_mass` via E²−p² loses precision for E ≫ m (use (E−p)(E+p)) |
| F34 | `timekeeping.rs:298` `AtomicInstant` serde | Deserialization bypasses nanosecond normalization, so `Ord` is wrong: {5 s, 4e9 ns} sorts before {8 s} |
| F35 | `optics.rs:40` `lines_to_spd` | `step_nm = 0` **panics** (integer overflow), violating the zero-panic rule; `fwhm = 0` gives NaN |
| F36 | `timekeeping.rs:391`, `:815` | Doc and test call −gΔh/c² a "blueshift"; it is a redshift. `second_order_doppler_shift` is −β²/2 (leading order) without saying so; −6.7% at β = 0.5 |

---

## Part 2 — Data tables

### D1. Provenance — most tables are not the edition they cite

| Table | Label says | Actually |
|---|---|---|
| `constants.rs` | CODATA 2022 | 10 of 13 non-exact values are **CODATA 2018**. `HBAR_MEV_S` is **2014** and differs by 8.4e-9 from the exact `HBAR_EV_S` two lines above it. `PROTON_G_FACTOR` is **2010** |
| `particle.rs` | PDG 2024 | m_d 4.67, m_s 93.4, m_c 1270, m_b 4180, m_H 125 250 are **PDG 2022/2023** (2024: 4.70, 93.5, 1273.0, 4183, 125 200) |
| `decay.rs` half-lives | NNDC / NUBASE2020 | **52 / 114** match NUBASE2020. 44 match only NUBASE2012/2016 or ENSDF, 4 only DDEP, and **13 match no evaluation checked** (incl. C-14, Cs-137, Th-232, U-234, Rn-222) |
| `reaction.rs` thermal | NNDC / Mughabghab 2018 | 13 of 16 absorption/scattering numbers are **Sears, Neutron News 3(3) 1992**: bound-atom scattering |
| `nucleus.rs` superallowed | Hardy & Towner 2020 | ft values: 8 of 9 identical to the **2009** survey. Average Ft 3072.27 ± 0.72 s is **2015** (2020: 3072.24 ± **1.85** s) |
| `nucleus.rs` moments | Stone 2005/2019 | 11 μ values are the **2014** table |
| `atomic.rs` ionization energies | "Z = 1–103 experimental, NIST ASD" | H is theory in NIST; Pa, Fm and Md are NIST interpolations; Es and No are Sugar 1974 estimates; Tc is a 1955 value |
| `atomic.rs` electron affinities | "NIST, Andersen 2004, Rienstra-Kiracofe 2002; Z > 104 Eliav/Kaldor/Borschevsky" | NIST ASD has no EA data (README says it does too). The superheavy attribution fits Rg, Nh and Ts only. Several leading digits match PubChem's periodic table exactly, with extra digits that match no source |
| `timekeeping.rs` Sr-87 / Yb-171 | CIPM 2021 | Sr-87 873.2 Hz is **CIPM 2015**; Yb-171 is **CIPM 2017**. CIPM 2025 values have been in force since 2026-03-27 |

### D2. Constants (`constants.rs`)

| Constant | Code | CODATA 2022 | Edition matched | Rel. diff |
|---|---|---|---|---|
| `ELECTRON_MASS_MEV` | 0.510 998 950 | 0.510 998 950 69 | 2018 | −1.35e-9 |
| `PROTON_MASS_MEV` | 938.272 088 16 | 938.272 089 43 | 2018 | −1.35e-9 |
| `NEUTRON_MASS_MEV` | 939.565 420 52 | 939.565 421 94 | 2018 | −1.51e-9 |
| `AMU_MEV` | 931.494 102 42 | 931.494 103 72 | 2018 | −1.40e-9 |
| `FINE_STRUCTURE` (1/α) | 137.035 999 084 | 137.035 999 177 | 2018 | α +6.8e-10 |
| `RYDBERG` | 10 973 731.568 160 | 10 973 731.568 157 | 2018 | +2.7e-13 |
| `BOHR_RADIUS` | 5.291 772 109 03e-11 | 5.291 772 105 44e-11 | 2018 | +6.8e-10 |
| `HBAR_MEV_S` | 6.582 119 514e-22 | 6.582 119 569…e-22 (exact) | **2014** | −8.4e-9 |
| `BOHR_MAGNETON_EV_T` | 5.788 381 8060e-5 | 5.788 381 7982e-5 | 2018 | +1.35e-9 |
| `CLASSICAL_ELECTRON_RADIUS_FM` | 2.817 940 3262 | 2.817 940 3205 | 2018 | +2.0e-9 |
| `PROTON_G_FACTOR` | 5.585 694 713 | 5.585 694 6893 | **2010** | +4.2e-9 |
| `ELECTRON_ANOMALOUS_MOMENT` | 1.159 652 181 28e-3 | 1.159 652 180 46e-3 | 2018 | +7.1e-10 |

The exact and defined constants (c, N_A, e, ħ in eV·s, h, k_B, g_n) and IERS `GM_EARTH`
and `EARTH_ROTATION_RAD_S` are correct. Two doc errors: `COULOMB_MEV_FM` says 1.4399764
(the value 1.4399645 is right), and `PROTON_G_FACTOR` says g_p = μ_p/μ_N (it is
2μ_p/μ_N). g_e in `atomic.rs` is labelled 2022 but is 2018.

Numerically all of this is ≤ 1e-8. It matters as a provenance claim, and because
constants are **re-hard-coded** outside this module:
- `197.3269804` in six places
- `13.6` three times
- `1239.842`
- `931.494`
- m_p, m_n and mₑ in `bridge.rs`

### D3. Nuclear structure (`nucleus.rs`)

| Sev | Item | Code | Reference | Note |
|---|---|---|---|---|
| **high** | F-19 electric quadrupole | −0.0942 b | **0**: ground state is 1/2⁺ (NUBASE2020) | −0.0942 b belongs to the 197 keV 5/2⁺ state. ⚠ LiveChart's F-19 row carries the same error |
| medium | AME2020 mass excess: **17 of 32 outside uncertainty** | Sn-120 −91105, Fe-56 −60601, Zr-90 −88768, Mo-98 −88113, U-235 40920.5, Ba-138 −88263, U-238 47308.9, S-32 −26016.16, P-31 −24440.99 | −91097.741, −60607.163, −88772.547, −88115.980, 40918.782, −88261.806, 47307.732, −26015.537, −24440.544 keV | up to 7.3 keV; `ame2020_fe56_mass_excess_negative` pins the wrong Fe-56 value |
| low | Ni-58 charge radius | 3.770 fm | 3.7757(20) | 2.9σ; the other 11 radii are exact |
| low | Bi-209 μ | 4.1106 μN | 4.092(2) (2019) | superseded |
| low | superallowed parents | Al-26, K-38 ground states; O-14 → N-14 g.s. | Al-26m (228 keV), K-38m (130 keV); O-14 feeds N-14 2.313 MeV 0⁺ | `(Z, A)` cannot express this |

### D4. Radioactive decay (`decay.rs`)

| Sev | Item | Code | NUBASE2020 | Note |
|---|---|---|---|---|
| **high** | Bi-212 primary mode | α | **β⁻ 64.06%**, α 35.94% | Th-232 chain follows the minor branch; Po-212 is unreachable |
| **high** | Ta-180m half-life and mode | 1.2e15 y, IT | **> 4.5e16 y**, no decay observed | an old *lower limit* (NUBASE2012) stored as a measurement, ≥ 37× too short |
| medium | Po-212 | 299 ns | 294.4(8) ns | +1.6%, 5.8σ |
| medium | Pb-214 | 26.8 min | 27.06(7) min | −0.96% |
| medium | Rb-87 | 4.923e10 y | 4.97(3)e10 y | −0.95% |
| medium | I-129 | 1.57e7 y | 1.614(12)e7 y | −2.7% |
| medium | Po-216, Po-214, Fe-55, Sr-90, Eu-155 | 0.145 s, 164.3 µs, 2.744 y, 28.79 y, 4.7611 y | 144.0 ms, 163.47 µs, 2.7562 y, 28.91 y, 4.742 y | 0.4–0.7%, 1.7–30σ |
| low | C-14 | 5730 y | 5700(30) y | old "Cambridge" value, pinned by a ±1 y test |
| low | year length | 365.25 d | 365.2422 d | +21 ppm on all 51 year-unit entries |
| low | Cu-64 primary mode | β⁺ | EC 44%, β⁻ 38.5%, β⁺ 17.5% | picks the smallest branch |
| low | Pa-234m and Ta-180m energies; 19 half-lives within 0.2% | | | older evaluations, truncation |

Chains: U-238 (14 steps) and U-235 (11) follow the main branches. Th-232 is correct
except at Bi-212. Two table entries (Pa-234, Tc-99m) are never used, because the
`(Z, A)` lookup takes the first match.

### D5. Reactions (`reaction.rs`)

| Sev | Item | Code | Reference (ENDF/B-VIII.0 unless noted) | Note |
|---|---|---|---|---|
| **high** | H-1 thermal scattering | 82.02 b | free atom **20.436 b** | 82.02 is Sears *bound* natural H: wrong quantity, ×4 |
| **high** | U-235 "capture" resonance integral | 275 b | capture **143.0 b** (fission 280.2 b) | it is the fission RI; pinned by a test |
| **high** | B-10 "capture" RI | 1722 b | (n,γ) **0.178 b**; (n,α) 1722.6 b | it is the (n,α) RI |
| **high** | Sn-120 RI | 0.133 b | **1.104 b** (Mughabghab 2003: 1.2(3)) | −88%; looks like thermal σ_γ |
| **high** | `u235_fission` Q for its own channel, n + U-235 → Ba-141 + Kr-92 + 3n | 200 MeV | **173.28 MeV** (AME2020) | 200 MeV is the total release including β decays: a different quantity |
| medium | `cno_cycle` Q | 25.03 MeV | Q = 4M(¹H)−M(⁴He) = **26.731 MeV** | 25.03 is Q minus mean ν losses; the pp preset uses the opposite convention |
| medium | Cd-113, C-12, U-235 thermal scattering | 12.1, 5.551, 15.04 b | 24.32, 4.748, 14.08 b | Sears coherent/bound values |
| medium | Pu-239 and Zr-90 RI | 200, 0.117 b | 180.1, 0.164 b | +11%, −29% |
| medium | Pu-239 chain yields A = 110, 134, 85 | 0.0040, 0.0708, 0.0054 | 0.00645, 0.07676, 0.00574 | −38%, −7.8%, −5.9% |
| medium | `NuclearReaction` data model | neutron stored as H-1 (`u235_fission`); products omit n, e⁺, ν | — | charge and baryon bookkeeping impossible |
| low | U-235 σ_a/σ_f, Pu-239 σ_a, Cd-113 σ_a, C-12 σ_a and RI, I-127/Cs-133/Mo-98 RI, U-235 yields A = 105/110/137/144, Pu-239 A = 106/147 | | | 0.3–22% |
| low | r-process "Os/Pt-195" | Z = 76 (Os-195, 6.5 min) | stable end product **Pt-195** (Z = 78) | |
| low | s-process | Fe-58 → "Fe-59→Co-59" | Fe-58(n,γ)Fe-59 step missing | |

⚠ **For regeneration:** the IAEA LiveChart fission-yield API serves **JEFF-3.1.1**, not
ENDF/B-VIII.0.

### D6. Timekeeping (`timekeeping.rs`, `bridge.rs`)

| Sev | Item | Code | Reference | Note |
|---|---|---|---|---|
| **high** | `fractional_stability`, "Allan deviation at 1 s" | H maser 1e-15; Sr/Yb 1e-18; Cs/Rb 1e-13 | H maser ~1.5e-13 (1e-15 only at ~10⁴ s); best Sr pair 4.8e-17 at 1 s; commercial Cs ≤ 1.2e-11, Rb < 2e-11 | 50–200× off, unsourced |
| medium | TCG/TCB conversions | `LG_RATE` and `LC_RATE` referenced only by tests; no functions; no L_B, TDB0 or T0 | IERS Conventions eqs. 10.1–10.5 | CHANGELOG ("with conversion functions") and roadmap mark this done |
| medium | GPS clock correction | 38.504 µs/day (R = 6371 km, no rotation) | 38.575 µs/day (IERS eq. 10.9, geoid potential) | docs quote 45.85 / 38.6, which the function does not return |
| low | Sr-87, Yb-171, Rb-87 frequencies | 873.2, 863.6, …904 312 Hz | CIPM 2025: 872.992, 863.632, …904 312 9 | Sr differs by 4.9e-16, resolvable in f64 |
| low | `quality_factor` | 1e10, 1e10, 1e9, 1e17, 1e17 | — | order-of-magnitude constants; the doc's Q = f·τ is never computed |
| low | labels | dTCG/dTT = 1+L_G; L_C "IAU 2006 B3"; mean radius "IERS 2010" | dTT/dTCG = 1−L_G; L_C from IERS Table 1.1; IERS defines no mean radius | |
| low | leap-second table | "as of Bulletin C 69" | Bulletin C 72 (2026-07-06): still 37 s; IANA file expires 2027-06-28 | no expiry handling |

Correct: all 28 leap-second entries, all 708 month lookups from 1972-01 to 2030-12,
JD/MJD epochs, Cs-133, H hyperfine, wavelengths, 32.184 s, 19 s, L_G, L_C, GM_earth
and Ω_earth.

### D7. Atomic data (`atomic.rs`)

| Sev | Item | Code | Reference | Note |
|---|---|---|---|---|
| **high** | Tc ionization energy | 7.28 eV | **7.11938(3)** (NIST; Mattolat 2010) | +0.161 eV, a 1955 value |
| **high** | Rg, Cn, Nh, Fl, Mc, Lv ionization energies | 9.79, 9.38, 5.85, 7.31, 6.92, 8.59 | 10.6, 11.97, 7.306, 8.539, 5.553, 6.881 (Eliav; Smits et al. 2023) | off by 0.8–2.6 eV; Nh–Lv look row-shuffled |
| **high** | Ds and Rg configurations | 6d⁹7s¹, 6d¹⁰7s¹ ("parallels Pt/Au") | **6d⁸7s², 6d⁹7s²** (relativistic predictions) | the "22 NIST exceptions" are NIST's 20 plus these 2 wrong ones |
| **high** | EA Pr, Th, Hf, Ga, Lu, Re, La, Tl | 0.962, 1.170, 0.017, 0.430, 0.340, 0.150, 0.470, 0.377 | 0.10923, 0.607690, 0.1780, 0.301166, 0.23882, 0.060396, 0.557546, 0.320053 | superseded by 2017–2020 measurements; 0.06–0.85 eV |
| **high** | EA U, Nd, Eu, Tb | 0.0 ("anion unstable") | 0.31497, 0.09748, 0.116, 0.13131 | measured bound anions reported as unbound |
| **high** | EA Fl, Mc | 0.905, 0.674 | Fl no bound anion; Mc 0.313–0.366 (theory) | |
| medium | Ds ionization energy | 8.7 | 9.56 (Smits 2023; predictions 9.5–11.2) | |
| medium | Cm, Es, Am, No ionization energies | 6.02196, 6.42, 5.9938, 6.65 | 5.992241, 6.36840, 5.97381, 6.62621 (NIST) | 0.02–0.05 eV; old estimates, Am digit typo |
| medium | EA Pt, Ni | 2.12810, 1.15616 | 2.12510, 1.15716 | digit typos |
| medium | EA Sr, Cr, Ce, In, Pb | 0.04816, 0.6660, 0.570, 0.404, 0.3643 | 0.05206, 0.67584, 0.600160, 0.38392, 0.356721 | superseded |
| medium | EA Lv; Gd, Db, Ds | 1.470; 0.0, 0.680, 0.0 | 0.776; 0.212, 1.189, 0.830 | theory, uncertain |
| low | Po, Pr ionization energies; Sg | 8.414, 5.473; 7.08 | 8.418070, [5.4702]; [7.8(5)] | old values |
| low | 20 ionization energies < 1 meV off | | NIST 2013 / CRC | Te, Ra, Lu and Sr carry digits no source has |
| low | EA Sc, Ti, V, Fe, Y, Zr, Nb, Mo, Ir (1–9 meV); Dy 0.0 (0.015 bound); Ag/Sb dropped digit; Og (0.076 calc.) | | | superseded |
| low | EA `0.0` semantics | means "unstable" (doc names alkaline earths, which have bound values) **and** silently "no data" (Pm, Sm, Ho, Er, Pa, Np, Am, Cm, Fm, Md, Rf, Bh, Hs, Mt) | | should be an explicit type |

Correct: 80 of 108 measured ionization energies within rounding and uncertainty;
116 of 118 configurations (all 20 NIST exceptions correct); 36 measured EAs; 14 unbound
anions correctly 0.0.

---

## Part 3 — Claims that are false

- **`docs/development/threat-model.md` numeric precision table:**
  - Rydberg "< 0.01%, exact for hydrogen-like" (F15: 0.054%).
  - Fine structure "~0.001%" (F16).
  - SEMF "~1–2%" (F20: RMS 10 MeV, biased).
  - Shell-corrected "~0.5–1% near magic numbers" (F8: worse everywhere).
  - Bateman "< 0.01%" (F10).
- **ADR-005** — the shell correction "improves accuracy for doubly-magic nuclei": it
  worsens them from +0.7 to +7.8 MeV.
- **`CHANGELOG.md` 1.2.0 / `roadmap.md`:**
  - "TCB, TCG … with conversion functions": the functions do not exist.
  - "Strutinsky shell correction" is a Gaussian heuristic.
  - "Breit interaction" is an unsourced formula.
- **`README.md` data sources:**
  - NIST ASD is listed for electron affinities; it has none.
  - "CODATA 2022", "PDG 2024", "NUBASE 2020" and "ENDF/B-VIII.0" labels (D1).
- **`CLAUDE.md` principles broken:**
  - "ALL physics values must use real data … No fake data, no stubs" — F8, F18,
    F22, F26, D6 `fractional_stability`/`quality_factor`, `optics.rs` Balmer
    intensities.
  - "Zero unwrap/panic" — F35.
  - "Never skip benchmarks" — 1.1/1.2 surface unbenchmarked.

---

## Part 4 — Implications for the Cyrius port

1. **Fix in Rust first, or port and fix together — a decision to make before
   planning.** The port will want the Rust crate as its oracle for golden values. As
   it stands the oracle is wrong in 36 formulas and ~100 table values. A Cyrius
   translation faithful to today's code would inherit all of it and bless it with
   parity tests.
2. **Regenerate every table from primary sources by script; do not hand-carry
   numbers.** The transcription error rate is too high: 17/32 AME rows, ~30/118 EA
   rows, 62/114 half-lives not the cited edition. Commit the generator, the source
   file edition and a checksum.
   - Sources: CODATA 2022 `allascii.txt`, AME2020 `mass_1.mas20`, NUBASE2020,
     NIST ASD CSV, ENDF/B-VIII.0, IERS `Leap_Second.dat` (with expiry), BIPM
     frequency API (CIPM 2025).
   - Two traps: LiveChart fission yields are JEFF-3.1.1, and LiveChart's F-19 row
     carries the D3 quadrupole error.
3. **Data model gaps to fix, not port:**
   - nuclide identity needs an isomer index (F11, D3, D4);
   - reactions need neutrons, leptons and neutrinos (D5);
   - decay data needs branching ratios, not a "primary mode" (D4 Bi-212, Cu-64);
   - EA needs an explicit "unbound" vs "no data" distinction and uncertainties (D7);
   - superallowed parents need isomers.
4. **Numerics to design in from the start:**
   - integer-nanosecond instant arithmetic, never through f64 (F13);
   - two-part Julian dates;
   - a stable Bateman formulation (divided differences / Φ-functions or matrix
     exponential) that handles equal λ (F10);
   - series branches for Klein–Nishina and the uniform form factor (F30, F31);
   - (E−p)(E+p) for invariant mass (F33);
   - validated inputs returning errors, not NaN, inf or panics (F32, F35);
   - reduced mass as an explicit spectroscopy input (F15).
5. **Replace heuristics with sourced physics, or delete them:**
   - shell correction (Myers–Swiatecki, or remove);
   - Lamb shift (tabulated F(Zα));
   - Einstein A (exact hydrogenic radial integrals / Gordon formula);
   - Breit (remove, or cite a real expression);
   - valence (derive from the configuration);
   - hyperfine scale;
   - soorat orbitals (`radial_wavefunction` × Y_lm);
   - clock stability and quality factor (cite datasheets/literature, or remove).
6. **Cyrius math surface:**
   - **No f64 cube root exists.** `pow(x, 1/3)` is wrong on 88% of perfect cubes and
     off by up to 137 ulps at the ends of the range. Filed in ganita as
     `ganita/docs/development/issues/2026-09-16-f64-cbrt-missing.md`, with a repro;
     measured on cyrius 6.6.4, identical under 6.6.2.
   - Cyrius exp/ln/pow will not match libm bit-for-bit, so golden tests need
     ulp-scale tolerances, not equality.
7. **Tests the port should carry over:**
   - reference values from primary data, not from the implementation;
   - property tests: normalization ∫R²r²dr = 1, atom conservation in Bateman, the
     optical theorem, the unscreened limit of Born → Rutherford, Z-scaling laws;
   - every finding here should have a test that **fails on the current code**.

## Reproduction

- **Formula findings:** a scratch Rust crate (path dependency on this repo, features
  `soorat-compat` and `optics`) prints each quantity. Independent references were
  computed alongside: Simpson integration, closed forms, and Python `decimal` at
  60–120 digits.
- **SEMF and spin-parity sweeps:** the crate's output for every (Z, A) was dumped and
  compared with AME2020 and NUBASE2020 in Python.
- **Data sweeps:** each table was parsed from the Rust source by regex and compared
  entry by entry with the primary files listed under Scope.

The harness and scripts live in the session scratch directory, so they are not
durable. Preserve them next to this report if the fixes should be re-verified
mechanically.

---

## Repairs — 2026-09-16 (v1.3.0)

Every finding above was repaired in the Rust crate the same day, so the Cyrius
port can take the Rust crate as its reference implementation. Each repair is
pinned by a test against the primary reference, not against the implementation.

| Findings | Repair | Pinned by |
|---|---|---|
| F1, F2 | Head-on turning point Z₁Z₂e²/T; total above θ_min = π(Z₁Z₂e²/2T)² cot²(θ/2) | `closest_approach_gold_alpha` (45.503 fm); `rutherford_total_matches_integral_of_differential` (1e-6) |
| F3 | Maximon (1968) full-range Bethe–Heitler | `pair_production_positive_and_continuous`, `pair_production_reference_values` |
| F4 | Exact R_nl for any n ≤ 60 (generalized Laguerre) | `radial_wavefunctions_normalized` (all n ≤ 10), `radial_r31_matches_textbook` |
| F5 | Exact hydrogenic dipole A from analytic radial integrals | `einstein_a_matches_nist_with_reduced_mass` (< 2e-4), `einstein_a_scales_as_z4` |
| F6, F18 | Uehling 4/15; one-loop Lamb shift with Bethe logarithms (`lamb_shift_nlj_ev`) | `vacuum_polarization_hydrogen_2s`, `lamb_shift_hydrogen_classic`, `lamb_shift_helium_ion` |
| F7 | Leading Breit–Pauli term α²Z³/4 E_h | `breit_leading_order_value`, `breit_scales_with_z3`, `breit_below_total_two_electron_energy_uranium` |
| F8, F20 | AME2020-fitted SEMF (RMS 3.31 MeV) + Myers–Swiatecki shell term (RMS 2.76 MeV) | `semf_close_to_ame2020`, `shell_correction_improves_doubly_magic`, `shell_term_sign_at_magic_and_midshell` |
| F9, F21 | Brennan–Bernstein coupling; empirical p/n level orders (odd-A 49%, odd-odd 30%, no half-integer spins) | `odd_odd_spin_is_integer`, `spin_parity_reference_nuclei` |
| F10 | Non-negative scaling-and-squaring matrix exponential | `bateman_equal_decay_constants`, `bateman_independent_of_time_unit`, `bateman_u238_series_trace_daughters` (80-digit references) |
| F11 | Ground state after IT; infinite half-life ends chain | `decay_chain_stops_after_isomeric_transition`, `bi212_follows_beta_branch` |
| F12, F26 | Valence from configuration; hyperfine factor I + ½ | `valence_from_configuration`, `hyperfine_hydrogen` |
| F13, F34 | Integer-nanosecond arithmetic; serde normalization | `atomic_instant_exact_at_modern_epoch`, `atomic_instant_serde_normalizes` |
| F14 | Closed-loop formula documented; new one-way `sagnac_one_way_ns` | `sagnac_one_way_equatorial_loop` (207.39 ns) |
| F15, F16 | `spectral_line_vacuum_nm`, `reduced_mass_factor`; CODATA Rydberg energy | `h_alpha_vacuum_with_reduced_mass` (656.4696 nm), `fine_structure_uses_codata_rydberg` |
| F17 | `hyperfine_splitting_spin_ev` with (I + ½) | `hyperfine_deuterium_includes_spin_factor` (327.384 MHz) |
| F19 | `atomic_mass_amu` includes Z mₑ − B_e; `nuclear_mass_amu` added | `atomic_mass_includes_electrons` |
| F22 | `corrected_ft_value` deprecated; per-transition `corrected_ft_seconds` | `superallowed_corrected_ft_per_transition` |
| F23, F24, F27, F28 | Mott from p·β; Mott factor, Paschen–Back and Born docs corrected; `born_screened_coulomb_with_masses` | `mott_electron_low_energy_uses_momentum` |
| F25 | Orbital slices use exact R_nl and \|Y_l0\|² | `orbital_*` soorat tests |
| F29 | `tai_minus_utc_seconds_mjd` covers 1961–1971 | `tai_minus_utc_rubber_second_era` |
| F30, F31 | Series branches (Klein–Nishina γ < 0.01, form factor qR < 0.5) | `klein_nishina_total_matches_integral_small_gamma` (1e-8), `form_factor_small_qr_no_cancellation` |
| F32, F33, F35, F36 | Input validation; precision-preserving relativity forms; `lines_to_spd` guards; redshift doc; `time_dilation_shift_exact` | `dirac_rejects_j_above_n`, `stark_and_lande_reject_impossible_states`, `gamma_symmetric_and_edge_cases`, `time_dilation_exact_vs_leading_order` |
| D1–D7 | Every table regenerated from its primary file: CODATA 2022, PDG 2024, AME2020, NUBASE2020, ENDF/B-VIII.0, Stone 2019/2021, Hardy & Towner 2020, NIST ASD, AHH99 + later EA measurements, Smits et al. 2023, CIPM 2025, IERS / USNO, Storey & Hummer 1995 | `ame2020_*`, `c14_half_life`, `u238_half_life`, `thermal_*`, `resonance_integral_*`, `fission_yields_match_endf_chain_yields`, `ionization_energy_reference_values`, `electron_affinity_reference_values`, `cipm_2025_frequencies`, `tcg_and_tdb_at_j2000`, `schwarzschild_gps_iers_reference` |

Clean after repair:
- `cargo fmt --check` and `cargo clippy --all-features --all-targets -D warnings`
- `cargo doc -D warnings`
- tests: 476 unit + 20 integration + 2 doc (all features); 464 + 20 + 2 (no default features)
- `cargo audit` and `cargo deny` (crossbeam-epoch 0.9.21)

Benchmarks: 19 suites, timings in `docs/guides/testing.md`. The one material
regression is `bateman_chain_3`, 88 ns → 1.0 µs, the price of the stable solver.

**Remaining model limits** (documented, not defects):
- the extreme single-particle shell model cannot describe deformed nuclei;
- the Lamb-shift Zα expansion is for low Z;
- pair production omits screening and the Coulomb correction;
- the Breit term is leading order only;
- tanmatra itself stores Ta-180m as non-decaying.

ganita shipped the f64 cube root filed above in **1.2.6**.
