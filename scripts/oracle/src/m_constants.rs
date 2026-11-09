//! `constants.tsv` — every public constant of `tanmatra::constants`, class A.
//!
//! These rows pin the generated `src/constants.cyr` (port plan §2 step 6) bit for bit.
//! Constants defined in other modules (`CNO_NEUTRINO_LOSS_MEV`, `LB_RATE`, ...) are
//! emitted with their own module.

use crate::tsv::{Class, Module, f, s};
use tanmatra::constants as k;

macro_rules! constants {
    ($m:ident; $($name:ident),* $(,)?) => {
        $( $m.row("const", Class::A, &[s(stringify!($name))], &[f(k::$name)]); )*
    };
}

pub fn emit() -> Module {
    let mut m = Module::new("constants");
    constants!(m;
        ELECTRON_MASS_MEV,
        PROTON_MASS_MEV,
        NEUTRON_MASS_MEV,
        AMU_MEV,
        FINE_STRUCTURE,
        HBAR_EV_S,
        RYDBERG,
        BOHR_RADIUS,
        AVOGADRO,
        ELEMENTARY_CHARGE,
        C,
        R0_FM,
        BOLTZMANN_EV,
        H_EV_S,
        HBAR_MEV_S,
        BOHR_MAGNETON_EV_T,
        COULOMB_MEV_FM,
        CLASSICAL_ELECTRON_RADIUS_FM,
        PROTON_G_FACTOR,
        ELECTRON_ANOMALOUS_MOMENT,
        STANDARD_GRAVITY,
        GM_EARTH,
        EARTH_ROTATION_RAD_S,
        HBAR_C_MEV_FM,
        HC_EV_NM,
        RYDBERG_EV,
        ELECTRON_MASS_U,
        PROTON_MASS_U,
        NEUTRON_MASS_U,
        HYDROGEN_ATOM_MASS_U,
        NUCLEAR_MAGNETON_EV_T,
        DEUTERON_G_FACTOR,
        ATOMIC_UNIT_TIME_S,
        COULOMB_K_SI,
    );
    m
}
