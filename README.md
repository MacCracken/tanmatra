# tanmatra

**tanmatra** (Sanskrit: subtle element) -- Atomic and subatomic physics library for the AGNOS project.

Standard Model particles, nuclear structure, radioactive decay, spectral lines, electron configurations, and nuclear reactions.

## Features

- **Standard Model**: All quarks, leptons, gauge bosons, and fundamental forces with real PDG masses
- **Nuclear Structure**: Bethe-Weizsacker binding energy, mass defect, nuclear radius, magic numbers
- **Radioactive Decay**: Alpha/beta/gamma decay, half-lives, activity, decay chains with 10 real isotopes
- **Atomic Structure**: Electron configurations (Aufbau + Cr/Cu exceptions), quantum number validation
- **Spectral Lines**: Rydberg formula, Balmer series for hydrogen
- **Ionization Energies**: NIST values for Z=1-36
- **Nuclear Reactions**: D-T/D-D fusion, p-p chain, U-235 fission, Coulomb barrier
- **Constants**: All CODATA 2022 fundamental constants

## Usage

```rust
use tanmatra::prelude::*;

// Iron-56: peak of binding energy curve
let bea = binding_energy_per_nucleon(IRON56.z, IRON56.a).unwrap();
assert!(bea > 8.5); // ~8.8 MeV/nucleon

// Hydrogen H-alpha spectral line
let h_alpha = balmer_series(3).unwrap();
assert!((h_alpha - 656.3).abs() < 0.5); // ~656.3 nm

// C-14 radioactive decay
let c14 = tanmatra::decay::carbon14();
let lambda = decay_constant(c14.half_life_s).unwrap();
// lambda ~ 3.836e-12 /s

// Electron configuration of chromium (exception)
let config = electron_configuration(24).unwrap();
// [Ar] 3d5 4s1
```

## Accuracy

- Physical constants: CODATA 2022 recommended values
- Ionization energies: NIST Atomic Spectra Database
- Quark/lepton masses: PDG 2024 review
- Nuclear half-lives: evaluated nuclear data (NNDC/ENSDF)
- Binding energies: Bethe-Weizsacker semi-empirical mass formula (approximate)

## Features

| Feature | Default | Description |
|---------|---------|-------------|
| `std` | yes | Standard library support |
| `logging` | no | `tracing` instrumentation |
| `full` | no | All features enabled |

## License

GPL-3.0-only
