//! # tanmatra
//!
//! Atomic and subatomic physics: Standard Model particles, nuclear structure,
//! radioactive decay, spectral lines, electron configurations, and nuclear reactions.
//!
//! All physical constants are from CODATA 2022. All ionization energies are from NIST.
//! The Bethe-Weizsacker semi-empirical mass formula is used for nuclear binding energies.
//! The Rydberg formula is used for spectral line wavelengths.
//!
//! ## Features
//!
//! - `std` (default) -- enables `std` support in serde and thiserror
//! - `logging` -- enables `tracing` instrumentation
//! - `full` -- enables all features
//!
//! ## Example
//!
//! ```
//! use tanmatra::nucleus::{binding_energy_per_nucleon, IRON56};
//! use tanmatra::atomic::balmer_series;
//!
//! // Iron-56: peak of binding energy curve
//! let bea = binding_energy_per_nucleon(IRON56.z, IRON56.a).unwrap();
//! assert!(bea > 8.5);
//!
//! // Hydrogen Balmer H-alpha line
//! let h_alpha = balmer_series(3).unwrap();
//! assert!((h_alpha - 656.3).abs() < 0.5);
//! ```

#![cfg_attr(not(feature = "std"), no_std)]
#![warn(missing_docs)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

extern crate alloc;

pub mod atomic;
pub mod constants;
pub mod decay;
pub mod error;
pub mod nucleus;
pub mod particle;
pub mod reaction;

/// Prelude: re-exports of the most commonly used types and functions.
pub mod prelude {
    pub use crate::atomic::{
        OrbitalType, QuantumNumbers, balmer_series, electron_configuration, ionization_energy_ev,
        max_electrons_in_shell, max_electrons_in_subshell, spectral_line_nm,
    };
    pub use crate::constants::*;
    pub use crate::decay::{
        DecayMode, KnownIsotope, activity_bq, alpha_decay, beta_minus_decay, beta_plus_decay,
        decay_chain, decay_constant, remaining_fraction,
    };
    pub use crate::error::TanmatraError;
    pub use crate::nucleus::{
        CARBON12, HELIUM4, HYDROGEN, IRON56, Nucleus, OXYGEN16, URANIUM235, URANIUM238,
        binding_energy_mev, binding_energy_per_nucleon, is_magic_number, mass_defect_mev,
        nuclear_radius_fm,
    };
    pub use crate::particle::{
        Boson, BosonProperties, ForceProperties, FundamentalForce, Lepton, LeptonProperties, Quark,
        QuarkProperties, boson_properties, force_properties, lepton_properties, quark_properties,
    };
    pub use crate::reaction::{
        NuclearReaction, coulomb_barrier_mev, dd_fusion, dt_fusion, is_exothermic, pp_chain_step1,
        q_value, u235_fission_approx,
    };
}

// Static assertions: all key types are Send + Sync.
#[cfg(test)]
const _: () = {
    #[allow(dead_code)]
    fn assert_send_sync<T: Send + Sync>() {}
    #[allow(dead_code)]
    fn _assertions() {
        assert_send_sync::<nucleus::Nucleus>();
        assert_send_sync::<particle::Quark>();
        assert_send_sync::<particle::Lepton>();
        assert_send_sync::<particle::Boson>();
        assert_send_sync::<particle::FundamentalForce>();
        assert_send_sync::<decay::DecayMode>();
        assert_send_sync::<decay::KnownIsotope>();
        assert_send_sync::<atomic::OrbitalType>();
        assert_send_sync::<atomic::QuantumNumbers>();
        assert_send_sync::<reaction::NuclearReaction>();
        assert_send_sync::<error::TanmatraError>();
    }
};
