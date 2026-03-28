//! Nuclear reactions: Q-values, fusion reactions, fission, Coulomb barrier.

extern crate alloc;
use alloc::vec;
use alloc::vec::Vec;

use crate::constants::COULOMB_CONST_MEV_FM;
use crate::nucleus::Nucleus;
use serde::{Deserialize, Serialize};

/// A nuclear reaction with reactants, products, and energy release.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NuclearReaction {
    /// Reactant nuclei.
    pub reactants: Vec<Nucleus>,
    /// Product nuclei.
    pub products: Vec<Nucleus>,
    /// Q-value in `MeV` (positive = exothermic).
    pub q_value_mev: f64,
}

/// Computes the Q-value of a nuclear reaction.
///
/// `Q = sum(reactant_masses) - sum(product_masses)`
#[must_use]
#[inline]
pub fn q_value(reactants_mass_mev: f64, products_mass_mev: f64) -> f64 {
    reactants_mass_mev - products_mass_mev
}

/// Returns true if the reaction is exothermic (Q > 0).
#[must_use]
#[inline]
pub fn is_exothermic(q: f64) -> bool {
    q > 0.0
}

/// D-T fusion: D + T -> He-4 + n, Q = 17.6.
///
/// The most energetically favorable fusion reaction, used in tokamak designs.
#[must_use]
pub fn dt_fusion() -> NuclearReaction {
    NuclearReaction {
        reactants: vec![
            Nucleus { z: 1, a: 2 }, // Deuterium
            Nucleus { z: 1, a: 3 }, // Tritium
        ],
        products: vec![
            Nucleus { z: 2, a: 4 }, // He-4
            Nucleus { z: 0, a: 1 }, // neutron (Z=0)
        ],
        q_value_mev: 17.6,
    }
}

/// D-D fusion: D + D -> He-3 + n, Q = 3.27.
///
/// One of two D-D branches; the other produces T + p.
#[must_use]
pub fn dd_fusion() -> NuclearReaction {
    NuclearReaction {
        reactants: vec![
            Nucleus { z: 1, a: 2 }, // Deuterium
            Nucleus { z: 1, a: 2 }, // Deuterium
        ],
        products: vec![
            Nucleus { z: 2, a: 3 }, // He-3
            Nucleus { z: 0, a: 1 }, // neutron
        ],
        q_value_mev: 3.27,
    }
}

/// Proton-proton chain step 1: p + p -> D + e+ + neutrino, Q = 1.44.
///
/// The dominant energy source in main-sequence stars like the Sun.
#[must_use]
pub fn pp_chain_step1() -> NuclearReaction {
    NuclearReaction {
        reactants: vec![
            Nucleus { z: 1, a: 1 }, // proton
            Nucleus { z: 1, a: 1 }, // proton
        ],
        products: vec![
            Nucleus { z: 1, a: 2 }, // Deuterium
                                    // positron and neutrino are not nuclei, but we account for them in Q
        ],
        q_value_mev: 1.44,
    }
}

/// Approximate U-235 fission: Q ~ 200.
///
/// Actual fission produces a distribution of products; this is the average
/// total energy release including kinetic energy and gamma rays.
#[must_use]
pub fn u235_fission_approx() -> NuclearReaction {
    NuclearReaction {
        reactants: vec![
            Nucleus { z: 92, a: 235 }, // U-235
            Nucleus { z: 0, a: 1 },    // neutron
        ],
        products: vec![
            // Approximate: Ba-141 + Kr-92 + 3n (one common channel)
            Nucleus { z: 56, a: 141 }, // Ba-141
            Nucleus { z: 36, a: 92 },  // Kr-92
            Nucleus { z: 0, a: 1 },    // neutron
            Nucleus { z: 0, a: 1 },    // neutron
            Nucleus { z: 0, a: 1 },    // neutron
        ],
        q_value_mev: 200.0,
    }
}

/// Coulomb barrier energy between two nuclei at distance r (in fm).
///
/// `V_c = k_e * Z1 * Z2 / r`
///
/// # Errors
///
/// Returns error if `r_fm` is zero, negative, NaN, or infinite.
pub fn coulomb_barrier_mev(
    z1: u16,
    z2: u16,
    r_fm: f64,
) -> Result<f64, crate::error::TanmatraError> {
    if r_fm <= 0.0 || r_fm.is_nan() || r_fm.is_infinite() {
        return Err(crate::error::TanmatraError::InvalidEnergy {
            reason: alloc::string::String::from("distance must be positive and finite"),
        });
    }
    Ok(COULOMB_CONST_MEV_FM * f64::from(z1) * f64::from(z2) / r_fm)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dt_fusion_q_value() {
        let r = dt_fusion();
        assert!((r.q_value_mev - 17.6).abs() < 0.1);
        assert!(is_exothermic(r.q_value_mev));
    }

    #[test]
    fn dd_fusion_q_value() {
        let r = dd_fusion();
        assert!((r.q_value_mev - 3.27).abs() < 0.1);
        assert!(is_exothermic(r.q_value_mev));
    }

    #[test]
    fn pp_chain_q_value() {
        let r = pp_chain_step1();
        assert!((r.q_value_mev - 1.44).abs() < 0.1);
    }

    #[test]
    fn u235_fission_q_value() {
        let r = u235_fission_approx();
        assert!((r.q_value_mev - 200.0).abs() < 1.0);
    }

    #[test]
    fn q_value_exothermic() {
        let q = q_value(1000.0, 990.0);
        assert!(is_exothermic(q));
        assert!((q - 10.0).abs() < f64::EPSILON);
    }

    #[test]
    fn q_value_endothermic() {
        let q = q_value(990.0, 1000.0);
        assert!(!is_exothermic(q));
    }

    #[test]
    fn coulomb_barrier_dt() {
        let r = 1.2 * (libm::cbrt(2.0) + libm::cbrt(3.0));
        let v = coulomb_barrier_mev(1, 1, r).unwrap();
        assert!(v > 0.3 && v < 0.6, "Coulomb barrier = {v} MeV");
    }

    #[test]
    fn coulomb_barrier_invalid_distance() {
        assert!(coulomb_barrier_mev(1, 1, 0.0).is_err());
        assert!(coulomb_barrier_mev(1, 1, -1.0).is_err());
    }

    #[test]
    fn conservation_dt_fusion() {
        let r = dt_fusion();
        let z_in: u16 = r.reactants.iter().map(|n| n.z).sum();
        let z_out: u16 = r.products.iter().map(|n| n.z).sum();
        let a_in: u16 = r.reactants.iter().map(|n| n.a).sum();
        let a_out: u16 = r.products.iter().map(|n| n.a).sum();
        assert_eq!(z_in, z_out, "charge conservation");
        assert_eq!(a_in, a_out, "baryon conservation");
    }

    #[test]
    fn conservation_u235_fission() {
        let r = u235_fission_approx();
        let z_in: u16 = r.reactants.iter().map(|n| n.z).sum();
        let z_out: u16 = r.products.iter().map(|n| n.z).sum();
        let a_in: u16 = r.reactants.iter().map(|n| n.a).sum();
        let a_out: u16 = r.products.iter().map(|n| n.a).sum();
        assert_eq!(z_in, z_out, "charge conservation");
        assert_eq!(a_in, a_out, "baryon conservation");
    }

    #[test]
    fn serde_roundtrip_reaction() {
        let r = dt_fusion();
        let json = serde_json::to_string(&r).unwrap();
        let back: NuclearReaction = serde_json::from_str(&json).unwrap();
        assert_eq!(r, back);
    }
}
