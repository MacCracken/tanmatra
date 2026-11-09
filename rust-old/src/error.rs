//! Error types for tanmatra.

/// Errors that can occur in tanmatra computations.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize, thiserror::Error)]
#[non_exhaustive]
pub enum TanmatraError {
    /// Invalid atomic number (Z must be >= 1).
    #[error("invalid atomic number Z={0}: must be >= 1")]
    InvalidAtomicNumber(u32),

    /// Invalid mass number (A must be >= Z).
    #[error("invalid mass number A={a} for Z={z}: A must be >= Z")]
    InvalidMassNumber {
        /// Atomic number.
        z: u32,
        /// Mass number.
        a: u32,
    },

    /// Invalid quantum numbers.
    #[error("invalid quantum numbers: {0}")]
    InvalidQuantumNumbers(alloc::string::String),

    /// Invalid half-life (must be positive and finite).
    #[error("invalid half-life: {0}")]
    InvalidHalfLife(alloc::string::String),

    /// Decay not possible for the given nucleus.
    #[error("decay not possible: {0}")]
    DecayNotPossible(alloc::string::String),

    /// Invalid reaction parameters.
    #[error("invalid reaction: {0}")]
    InvalidReaction(alloc::string::String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn serde_roundtrip_error() {
        let err = TanmatraError::InvalidAtomicNumber(42);
        let json = serde_json::to_string(&err).unwrap();
        let back: TanmatraError = serde_json::from_str(&json).unwrap();
        assert_eq!(err, back);
    }

    #[test]
    fn serde_roundtrip_error_with_string() {
        let err = TanmatraError::InvalidQuantumNumbers(alloc::string::String::from("test"));
        let json = serde_json::to_string(&err).unwrap();
        let back: TanmatraError = serde_json::from_str(&json).unwrap();
        assert_eq!(err, back);
    }
}
