//! # tanmatra — Atomic and Subatomic Physics
//!
//! **tanmatra** (Sanskrit: तन्मात्र — subtle element) provides atomic and
//! subatomic physics computations: Standard Model particles, nuclear structure,
//! radioactive decay, spectral lines, and nuclear reactions.
//!
//! ## Architecture
//!
//! ```text
//! Fundamental ─── particle.rs ─── Standard Model (quarks, leptons, bosons, forces)
//!      |
//!      v
//! Nuclear ──────── nucleus.rs ─── Bethe-Weizsacker binding energy, nuclear radii
//!      |            decay.rs ──── Radioactive decay, half-lives, decay chains
//!      |            reaction.rs ─ Nuclear reactions, Q-values, Coulomb barriers
//!      v
//! Atomic ────────── atomic.rs ─── Electron configurations, spectral lines, ionization
//! ```
//!
//! ## Quick Start
//!
//! ```rust
//! use tanmatra::prelude::*;
//!
//! // Iron-56 binding energy per nucleon
//! let fe56 = Nucleus::iron_56();
//! let bea = fe56.binding_energy_per_nucleon();
//! assert!(bea > 8.4 && bea < 9.2); // ~8.8 MeV
//!
//! // H-alpha spectral line
//! let h_alpha = spectral_line_nm(1, 2, 3).unwrap();
//! assert!((h_alpha - 656.3).abs() < 1.0);
//!
//! // Electron configuration of iron
//! let config = electron_configuration(26).unwrap();
//! let short = format_configuration_short(&config, 26);
//! assert_eq!(short, "[Ar] 4s2 3d6");
//! ```
//!
//! ## Feature Flags
//!
//! | Feature | Default | Description |
//! |---------|---------|-------------|
//! | `std` | Yes | Standard library support. Disable for `no_std` + `alloc` |
//! | `logging` | No | Structured tracing via the `tracing` crate |
//! | `full` | No | Enables all optional features |
//!
//! ## Data Sources
//!
//! - **CODATA 2022**: Fundamental physical constants (NIST)
//! - **PDG 2024**: Particle masses (Particle Data Group)
//! - **NNDC/NUBASE**: Nuclear half-lives (National Nuclear Data Center)
//! - **NIST ASD**: Ionization energies (Atomic Spectra Database)

#![cfg_attr(not(feature = "std"), no_std)]
#![forbid(unsafe_code)]
#![warn(missing_docs, clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::similar_names)]
#![allow(clippy::unreadable_literal)]
#![allow(clippy::doc_markdown)]
#![allow(clippy::cast_lossless)]
#![allow(clippy::match_same_arms)]

extern crate alloc;

pub mod atomic;
pub mod constants;
pub mod decay;
pub mod error;
pub mod nucleus;
pub mod particle;
pub mod reaction;

/// Prelude module — import everything commonly needed.
pub mod prelude {
    pub use crate::atomic::{
        OrbitalFilling, OrbitalType, QuantumNumbers, TransitionType, check_selection_rules,
        check_selection_rules_full, einstein_a_coefficient, einstein_b_coefficient,
        electron_affinity_ev, electron_configuration, format_configuration,
        format_configuration_short, hydrogen_level_energy_ev, ionization_energy_ev, lande_g_factor,
        radial_probability_density, radial_wavefunction, spectral_line_fine_nm, spectral_line_nm,
        stark_shift_hydrogen_ev, zeeman_splitting_ev,
    };
    pub use crate::constants::*;
    pub use crate::decay::{
        DecayMode, Isotope, activity_bq, alpha_decay, bateman_chain, beta_minus_decay,
        beta_plus_decay, decay_chain, decay_constant, known_isotopes, remaining_fraction,
    };
    pub use crate::error::TanmatraError;
    pub use crate::nucleus::{
        Nucleus, ShellLevel, ground_state_spin_parity, is_magic_number, next_shell_closure,
        shell_closure_below, shell_model_levels, shell_occupation,
    };
    pub use crate::particle::{Boson, FundamentalForce, Lepton, Quark};
    pub use crate::reaction::{
        FissionYield, NuclearReaction, NucleosynthesisPathway, NucleosynthesisProcess,
        NucleosynthesisStep, ThermalCrossSection, average_lethargy_gain,
        breit_wigner_cross_section, cno_cycle, collisions_to_thermalize, coulomb_barrier,
        dd_fusion_he3, dd_fusion_t, dt_fusion, geometric_cross_section_barns,
        max_energy_loss_fraction, moderating_ratio, pp_chain_step1, pu239_fission_yields, q_value,
        r_process_main, s_process_main, thermal_neutron_cross_sections, triple_alpha, u235_fission,
        u235_fission_yields,
    };
}

// Static assertions: all key types are Send + Sync.
#[allow(dead_code)]
trait AssertSendSync: Send + Sync {}
impl AssertSendSync for nucleus::Nucleus {}
impl AssertSendSync for particle::Quark {}
impl AssertSendSync for particle::Lepton {}
impl AssertSendSync for particle::Boson {}
impl AssertSendSync for particle::FundamentalForce {}
impl AssertSendSync for decay::DecayMode {}
impl AssertSendSync for atomic::OrbitalType {}
impl AssertSendSync for atomic::QuantumNumbers {}
impl AssertSendSync for atomic::OrbitalFilling {}
impl AssertSendSync for error::TanmatraError {}
