//! `constants.tsv` — every public constant of `tanmatra::constants`, class A.
//!
//! These rows pin the generated `src/constants.cyr` bit for bit. The list is shared with
//! that generator (`t_constants::all`), so the two can't drift apart. Constants defined in
//! other modules (`CNO_NEUTRINO_LOSS_MEV`, `LB_RATE`, ...) are emitted with their own
//! module.

use crate::t_constants;
use crate::tsv::{Class, Module, f, s};

pub fn emit() -> Module {
    let mut m = Module::new("constants");
    for &(name, value) in t_constants::all() {
        m.row("const", Class::A, &[s(name)], &[f(value)]);
    }
    m
}
