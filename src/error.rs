//! Error types for brahmanda.

use thiserror::Error;

/// Errors that can occur in large-scale structure computations.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum BrahmandaError {
    /// Invalid galaxy parameters (mass, radius, etc.).
    #[error("invalid galaxy parameter: {0}")]
    InvalidGalaxy(String),

    /// Invalid cosmological structure parameter.
    #[error("invalid structure parameter: {0}")]
    InvalidStructure(String),

    /// Numerical computation failed to converge.
    #[error("convergence failure after {iterations} iterations: {detail}")]
    ConvergenceFailed { iterations: usize, detail: String },

    /// Generic computation error (overflow, division by zero, etc.).
    #[error("computation error: {0}")]
    Computation(String),

    /// Non-finite input (NaN or Infinity).
    #[error("non-finite input in {context}: {value}")]
    NonFinite { context: &'static str, value: f64 },
}

/// Validate that a single f64 input is finite.
#[inline]
pub fn require_finite(value: f64, context: &'static str) -> Result<(), BrahmandaError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(BrahmandaError::NonFinite { context, value })
    }
}

/// Validate that all f64 inputs are finite.
#[inline]
pub fn require_all_finite(values: &[f64], context: &'static str) -> Result<(), BrahmandaError> {
    for &v in values {
        require_finite(v, context)?;
    }
    Ok(())
}

/// Validate that a computed result is finite.
#[inline]
pub fn ensure_finite(value: f64, context: &'static str) -> Result<f64, BrahmandaError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(BrahmandaError::Computation(format!(
            "{context}: result is {value}"
        )))
    }
}
