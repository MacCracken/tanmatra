//! `src/constants.cyr` — every public constant of `tanmatra::constants`, as a `TM_` const.
//!
//! The values come from the oracle's API. The comments are the Rust doc comments, which
//! carry the CODATA 2022 citations. `tests/golden/constants.tsv` pins every bit.

use crate::cyr::{doc_comment, f64_lit, header};
use tanmatra::constants as k;

macro_rules! constants {
    ($($name:ident),* $(,)?) => {
        &[ $( (stringify!($name), k::$name) ),* ]
    };
}

const ALL: &[(&str, f64)] = constants!(
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

/// The (name, value) pairs, for the golden emitter too.
pub fn all() -> &'static [(&'static str, f64)] {
    ALL
}

pub fn emit(src: &str, provenance: &str) -> Result<String, String> {
    let mut out = header(
        "physical constants (CODATA 2022 unless noted)",
        provenance,
        &[
            "Each value is the shortest decimal that round-trips to the oracle's f64, written out\n\
           positionally: the Cyrius lexer has no exponent notation, and decimal literals are\n\
           correctly rounded (6.6.10+). tests/golden/constants.tsv pins every bit.",
        ],
    );
    for &(name, value) in ALL {
        out.push_str(&doc_comment(src, &format!("pub const {name}:"))?);
        out.push_str(&format!("const TM_{name} = {};\n\n", f64_lit(value)?));
    }
    Ok(out.trim_end().to_string() + "\n")
}
