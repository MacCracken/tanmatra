//! Error types for tanmatra.

extern crate alloc;

/// Errors that can occur in tanmatra computations.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize, thiserror::Error)]
#[non_exhaustive]
pub enum TanmatraError {
    /// Invalid element: atomic number out of range (1-118).
    #[error("invalid element: Z={z} is out of range 1-118")]
    InvalidElement {
        /// The invalid atomic number.
        z: u16,
    },

    /// Invalid isotope: mass number inconsistent with atomic number.
    #[error("invalid isotope: Z={z}, A={a} (A must be >= Z and > 0)")]
    InvalidIsotope {
        /// Atomic number.
        z: u16,
        /// Mass number.
        a: u16,
    },

    /// Invalid particle specification.
    #[error("invalid particle: {reason}")]
    InvalidParticle {
        /// Description of why the particle is invalid.
        reason: alloc::string::String,
    },

    /// Decay operation failed.
    #[error("decay failed: {reason}")]
    DecayFailed {
        /// Description of why the decay failed.
        reason: alloc::string::String,
    },

    /// Invalid energy value.
    #[error("invalid energy: {reason}")]
    InvalidEnergy {
        /// Description of the energy error.
        reason: alloc::string::String,
    },

    /// General computation error.
    #[error("computation error: {reason}")]
    ComputationError {
        /// Description of the computation error.
        reason: alloc::string::String,
    },
}
