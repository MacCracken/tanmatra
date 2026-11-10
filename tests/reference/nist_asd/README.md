# NIST ASD energy levels: He II and Li III

These are the reference data of `tests/hydrogenic_nist.tcyr`, which checks tanmatra's hydrogen-like functions at Z = 2
and 3. The two `.tsv` files are kept byte for byte as NIST served them. `scripts/reference/nist_levels.py` turns them
into `tests/nist_levels_data.cyr`. Don't edit either file by hand: download a new snapshot and regenerate.

## Source

Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2024). *NIST Atomic Spectra Database* (ver. 5.12),
https://physics.nist.gov/asd. Retrieved 2026-11-09.

| File | Spectrum | Isotope (NIST's note) | Levels | Ionization limit, cm⁻¹ | SHA-256 |
|---|---|---|---|---|---|
| `he_ii_levels.tsv` | He II | ⁴He | 148, n = 1–20 | 438908.878840(8) | `b5af0e773fc79e6af11c76046e85c3419ad0b149dd3f85c8b6d933610a1e870a` |
| `li_iii_levels.tsv` | Li III | ⁷Li | 148, n = 1–20 | 987661.0139(9) | `cd1eec6c1d617b114fcbe50e97ac9b42ac54242b4134b5ec0d9191a05b0554a2` |

Every (n, l, j) level is listed for n ≤ 11. For n = 12–20, NIST lists only nS½, nP½ and the level with the highest l.
Every level is a calculated value, which ASD marks with parentheses. They come from two sources:
- **n = 1, n = 2 and the ionization limit:** V. A. Yerokhin and V. M. Shabaev, "Lamb Shift of n = 1 and n = 2 States
  of Hydrogen-like Atoms, 1 ≤ Z ≤ 110", *J. Phys. Chem. Ref. Data* **44**, 033103 (2015),
  doi:10.1063/1.4927487 (ASD reference L19200).
- **n ≥ 3:** G. W. Erickson, "Energy Levels of One-Electron Atoms", *J. Phys. Chem. Ref. Data* **6**, 831–869 (1977),
  doi:10.1063/1.555557 (ASD reference L3620, comment c107). NIST's note, verbatim: "Erickson's theoretical values of
  ionization energies were scaled to the CODATA-2018 values of the fundamental constants and then subtracted from the
  ground-state ionization energy given by Yerokhin & Shabaev 2015. Additional corrections were applied to the nS
  levels, as described by Kramida et al. 2008."

Levels are in cm⁻¹ above the ground level 1s ²S½. NIST states the uncertainty of the n = 2 levels as an excitation
energy from the ground level, and of every other level as a separation from the ionization limit.

## Query

Both files came from ASD's energy-level form with tab-delimited output (`format=3`), every column on, and the
bibliography (`biblio=on`). The He II URL was:

```
https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&spectrum=He+II&units=0&format=3&output=0&page_size=15&multiplet_ordered=0&conf_out=on&term_out=on&level_out=on&unc_out=1&j_out=on&biblio=on&temp=&submit=Retrieve+Data
```

The Li III URL is the same with `spectrum=Li+III`.
