# tanmatra

**tanmatra** (Sanskrit: तन्मात्र — subtle element) — Atomic and subatomic physics for [AGNOS](https://github.com/MacCracken/agnosticos).

[![CI](https://github.com/MacCracken/tanmatra/actions/workflows/ci.yml/badge.svg)](https://github.com/MacCracken/tanmatra/actions/workflows/ci.yml)
[![license](https://img.shields.io/badge/license-GPL--3.0--only-blue.svg)](LICENSE)

tanmatra **2.0.2** is a Cyrius library: a parity port of the Rust crate **1.3.0** (the crates.io and docs.rs
release), which stays frozen in [`rust-old/`](rust-old/) as the reference implementation. Every function is
checked against it value by value:
- **191,583 golden rows**, all within their tolerance class and almost all bit-identical, on x86_64 Linux,
  Windows PE and aarch64 ([port parity](docs/development/port-parity.md));
- **the Rust tests**: every one of the 498 is ported or accounted for in the
  [Rust-test ledger](docs/development/rust-test-ledger.md).

Beyond Rust, the hydrogen-like functions are also checked against the NIST He II and Li III energy levels at Z = 2
and 3 ([reference validation](docs/development/reference-validation.md)).

The [port plan](docs/development/tanmatra-port-plan.md) records the decisions, and [`state.md`](docs/development/state.md)
the current state.

## Modules

| Module (`src/`) | Description |
|--------|-------------|
| `constants` | CODATA 2022 fundamental physical constants (`TM_*`) |
| `particle` | Standard Model: quarks, leptons, bosons, forces (PDG 2024) |
| `nucleus`, `shell_model` | AME2020-fitted binding energy + Myers–Swiatecki shell term, shell model, AME2020 masses, radii, moments, superallowed Ft |
| `decay` | Radioactive decay, 114 isotopes (NUBASE2020), stable Bateman solver, decay chains |
| `atomic_config`, `atomic_hydrogenic`, `atomic_qed` | Electron configs, spectral lines, ionization/affinity, exact hydrogenic wavefunctions and A-values, Zeeman/Stark, Lamb shift, Dirac |
| `reaction` | AME2020 Q-values, ENDF/B-VIII.0 thermal cross sections, resonance integrals, fission yields, moderation, nucleosynthesis |
| `relativity` | Four-momentum, Lorentz factor, velocity addition, de Broglie |
| `scattering` | Rutherford/Mott/Born cross-sections, form factors, Klein–Nishina, Bethe–Heitler pair production |
| `timekeeping` | CIPM 2025 frequency standards, TAI/UTC/TT/GPS/TCG/TCB/TDB, leap seconds, relativistic clock corrections |
| `optics` | prakash spectral power distributions and colours (the optional `optics` profile) |
| `bridge` | The chrono bridge: ΔAT for a UTC year and month, TAI ↔ UTC seconds since 1958 (a `lib/chrono.cyr` datetime supplies them, from 1970 on). The other bridges wait on their consumers |
| `error` | `TM_ERR_*` codes for the 6 `TanmatraError` variants |

The tables (`*_data.cyr`, `constants.cyr`) are generated from the oracle by [`scripts/oracle/`](scripts/oracle/),
never typed by hand.

## Using it

Declare it in your `cyrius.cyml` and vendor it with `cyrius deps`:

```toml
[deps.tanmatra]
git = "https://github.com/MacCracken/tanmatra.git"
tag = "2.0.2"
modules = ["dist/tanmatra.cyr"]
```

```cyrius
include "lib/tanmatra.cyr"

fn main(): i64 {
    alloc_init();

    # Iron-56 binding energy per nucleon (~8.8 MeV)
    var bea: f64 = nucleus_binding_energy_per_nucleon(nucleus_iron_56());

    # H-alpha vacuum wavelength with the proton reduced mass (656.470 nm)
    var t, v = spectral_line_vacuum_nm(1, 1, 2, 3);
    var h_alpha: f64 = v;

    # Electron configuration of iron: [Ar] 4s2 3d6
    var cfg: TmOrbitalFilling[20];
    var ct, n = electron_configuration(26, &cfg, 20);
    var short = format_configuration_short(&cfg, n, 26);

    # A proton at 500 MeV/c
    var p = TmFourMomentum_from_mass_and_momentum(TM_PROTON_MASS_MEV, 500.0);
    var gamma: f64 = p.gamma();
    return 0;
}
```

**The API follows Rust's** (port plan §5 has the full mapping):
- **Names.**
  - Free functions keep their Rust names.
  - An enum's or a packed value's methods are `snake(type)_method(v)`, as in `nucleus_binding_energy(n)`.
  - Structs are `Tm<Name>` with an `impl`, as in `TmFourMomentum`, `TmAtomicInstant` and `TmNuclearReaction`.
  - Globals and enum members carry `TM_`.
- **Results and options.**
  - A fallible function returns a value-form `Result` (`var t, v = f(..); is_ok(t)`) whose Err payload is a `TM_ERR_*`
    code.
  - An optional value is a value-form `Option`.
- **Collections.** Rust's `Vec` returns fill a caller array and return the count.

The `optics` profile adds `src/optics.cyr` and needs prakash 2.15.2. Declare `[deps.prakash]` (tag 2.15.2) beside
tanmatra, list `"sakshi"` in `[deps] stdlib`, and use `modules = ["dist/tanmatra-optics.cyr"]`.

## Building

Cyrius (toolchain pinned in `cyrius.cyml`):

```bash
cyrius deps --features optics             # vendor the stdlib and prakash into lib/
cyrius build src/main.cyr build/tanmatra  # compile
cyrius test                               # tests/*.tcyr: the parity harnesses, the ported Rust tests and the NIST checks
cyrius bench                              # tests/*.bcyr
cyrius audit                              # fmt + lint + test + bench
cyrius distlib --all                      # dist/tanmatra.cyr and dist/tanmatra-optics.cyr
```

The oracle and its checks:

```bash
cargo test --manifest-path rust-old/Cargo.toml --all-features                       # the Rust 1.3.0 oracle, read-only
cargo run --release --locked --manifest-path scripts/oracle/Cargo.toml -- --check   # golden corpus and tables current
python3 -I scripts/surface_parity.py .                                              # every Rust item has its Cyrius name
python3 -I scripts/rust_test_ledger.py . --check                                    # every Rust test is accounted for
```

## Data Sources

| Source | Used For |
|--------|----------|
| [CODATA 2022](https://physics.nist.gov/cuu/Constants/) | Fundamental constants |
| [PDG 2024](https://pdg.lbl.gov/) | Particle masses, decay widths |
| [AME2020 / NUBASE2020](https://www-nds.iaea.org/amdc/) | Atomic masses; half-lives, decay modes, J^π (114 isotopes) |
| [NIST ASD](https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html) | Ionization energies (Z ≤ 108), hydrogen A-values |
| Andersen, Haugen & Hotop (1999) + later measurements; Smits et al. (2023) | Electron affinities; superheavy ionization energies |
| [ENDF/B-VIII.0](https://www.nndc.bnl.gov/endf/) | Thermal cross sections, resonance integrals, fission yields |
| Angeli & Marinova (2013); Stone INDC(NDS)-0794/0833 | Charge radii; magnetic dipole and quadrupole moments |
| Hardy & Towner (2020) | Superallowed ft and Ft values |
| CIPM 2025, IERS Conventions 2010, IERS Bulletin C | Frequency standards, time scales, leap seconds |
| Storey & Hummer (1995) | Case B Balmer intensities |

## Features

| Feature | Default | Description |
|---------|---------|-------------|
| `optics` | No | prakash spectral power distributions and colours: `[lib.optics]` and `dist/tanmatra-optics.cyr` |

The Rust crate's `std`, `logging` and `full` features have no Cyrius counterpart. The bridges and the soorat
integration wait on their consumers (port plan §6).

## License

GPL-3.0-only. See [LICENSE](LICENSE).
