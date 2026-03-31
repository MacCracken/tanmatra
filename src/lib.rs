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
/// Cross-crate bridges — primitive-value conversions from other AGNOS science crates.
pub mod bridge;
pub mod constants;
pub mod decay;
pub mod error;
/// Integration APIs for downstream consumers (soorat rendering).
pub mod integration;
pub mod nucleus;
pub mod particle;
pub mod reaction;
pub mod relativity;
pub mod scattering;
/// Frequency standards, atomic time scales, and relativistic clock corrections.
pub mod timekeeping;

/// Optics integration with the prakash crate.
///
/// Requires the `optics` feature flag: `tanmatra = { features = ["optics"] }`.
#[cfg(feature = "optics")]
pub mod optics;

/// Prelude module — import everything commonly needed.
pub mod prelude {
    pub use crate::atomic::{
        OrbitalFilling, OrbitalType, QuantumNumbers, TransitionType, anomalous_zeeman_splitting_ev,
        balmer_series, bound_electron_g_factor, brackett_series, breit_interaction_ev,
        check_selection_rules, check_selection_rules_full, dirac_binding_energy_ev,
        dirac_energy_mev, einstein_a_coefficient, einstein_b_coefficient, electron_affinity_ev,
        electron_configuration, electron_g_factor, format_configuration,
        format_configuration_short, hydrogen_level_energy_ev, hyperfine_splitting_ev,
        ionization_energy_ev, lamb_shift_ev, lande_g_factor, lyman_series, paschen_series,
        pfund_series, radial_probability_density, radial_wavefunction, relativistic_correction_ev,
        spectral_line_fine_nm, spectral_line_nm, stark_shift_hydrogen_ev, vacuum_polarization_ev,
        zeeman_splitting_ev,
    };
    pub use crate::bridge::{SimulationClock, TimeContext};
    pub use crate::constants::*;
    pub use crate::decay::{
        DecayMode, Isotope, activity_bq, alpha_decay, bateman_chain, beta_minus_decay,
        beta_plus_decay, decay_chain, decay_constant, known_isotopes, remaining_fraction,
    };
    pub use crate::error::TanmatraError;
    pub use crate::nucleus::{
        NuclearMoments, Nucleus, ShellLevel, SuperallowedDecay, corrected_ft_value,
        ground_state_spin_parity, is_magic_number, next_shell_closure, shell_closure_below,
        shell_model_levels, shell_occupation, superallowed_ft_values,
    };
    pub use crate::particle::{Boson, FundamentalForce, Lepton, Quark};
    pub use crate::reaction::{
        FissionYield, NuclearReaction, NucleosynthesisPathway, NucleosynthesisProcess,
        NucleosynthesisStep, ThermalCrossSection, average_lethargy_gain,
        breit_wigner_cross_section, cno_cycle, collisions_to_thermalize, coulomb_barrier,
        dd_fusion_he3, dd_fusion_t, dt_fusion, geometric_cross_section_barns,
        max_energy_loss_fraction, moderating_ratio, pp_chain_step1, pu239_fission_yields, q_value,
        r_process_main, resonance_integral_barns, s_process_main, thermal_neutron_cross_sections,
        triple_alpha, u235_fission, u235_fission_yields,
    };
    pub use crate::relativity::{
        FourMomentum, de_broglie_wavelength_fm, gamma_to_beta, invariant_mass_two_body,
        lorentz_gamma, relativistic_energy, relativistic_momentum, velocity_addition,
        velocity_to_beta,
    };
    pub use crate::scattering::{
        born_screened_coulomb, compton_energy_ratio, distance_of_closest_approach,
        klein_nishina_differential, klein_nishina_total, legendre_polynomial,
        mott_correction_factor, mott_electron_differential, mott_electron_with_form_factor,
        nuclear_form_factor_uniform, pair_production_cross_section, partial_wave_cross_section,
        partial_wave_differential, rutherford_differential, rutherford_total_above_angle,
        sommerfeld_parameter, thomas_fermi_screening_fm,
    };
    pub use crate::timekeeping::{
        AtomicInstant, FrequencyStandard, TimeScale, all_frequency_standards, gps_to_tai,
        gravitational_redshift, leap_seconds_at, sagnac_correction_ns,
        schwarzschild_clock_correction_us_per_day, second_order_doppler_shift, tai_to_gps,
        tai_to_tt, tai_to_utc_offset, tt_to_tai, utc_to_tai_offset,
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
impl AssertSendSync for relativity::FourMomentum {}
impl AssertSendSync for atomic::TransitionType {}
impl AssertSendSync for nucleus::ShellLevel {}
impl AssertSendSync for nucleus::NuclearMoments {}
impl AssertSendSync for nucleus::SuperallowedDecay {}
impl AssertSendSync for timekeeping::FrequencyStandard {}
impl AssertSendSync for timekeeping::TimeScale {}
impl AssertSendSync for timekeeping::AtomicInstant {}
impl AssertSendSync for bridge::SimulationClock {}
impl AssertSendSync for bridge::TimeContext {}
