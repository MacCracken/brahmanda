//! Galaxy morphology — Hubble classification, structural parameters.
//!
//! Classification of galaxies by type and measurement of their
//! structural properties (effective radius, Sersic profile, mass-to-light).

use serde::{Deserialize, Serialize};

use crate::error::{ensure_finite, require_finite, BrahmandaError};

/// Hubble classification of galaxy morphology.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum HubbleType {
    /// Elliptical (E0–E7).
    Elliptical,
    /// Lenticular (S0).
    Lenticular,
    /// Spiral (Sa, Sb, Sc, Sd).
    Spiral,
    /// Barred spiral (SBa, SBb, SBc).
    BarredSpiral,
    /// Irregular.
    Irregular,
}

/// Structural parameters of a galaxy.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[must_use]
pub struct GalaxyProperties {
    /// Hubble morphological type.
    pub hubble_type: HubbleType,
    /// Total stellar mass (solar masses).
    pub stellar_mass_msun: f64,
    /// Effective (half-light) radius (kpc).
    pub effective_radius_kpc: f64,
    /// Sersic index n (1=exponential disk, 4=de Vaucouleurs elliptical).
    pub sersic_index: f64,
    /// Velocity dispersion (km/s) — for ellipticals / bulges.
    pub velocity_dispersion_km_s: f64,
}

/// Sersic surface brightness profile: I(R) = I_e exp{-b_n [(R/R_e)^(1/n) - 1]}.
///
/// Returns the intensity ratio I(R)/I_e.
///
/// # Arguments
/// * `r` — Galactocentric radius.
/// * `r_e` — Effective (half-light) radius (same units as r).
/// * `n` — Sersic index (must be > 0).
///
/// ```
/// use brahmanda::morphology::sersic_profile;
///
/// // At R = R_e, the profile returns ~1.0
/// let i = sersic_profile(1.0, 1.0, 4.0).unwrap();
/// assert!((i - 1.0).abs() < 0.01);
///
/// // Intensity decreases outward
/// let i_outer = sersic_profile(3.0, 1.0, 4.0).unwrap();
/// assert!(i_outer < i);
/// ```
#[inline]
pub fn sersic_profile(r: f64, r_e: f64, n: f64) -> Result<f64, BrahmandaError> {
    require_finite(r, "sersic_profile")?;
    require_finite(r_e, "sersic_profile")?;
    require_finite(n, "sersic_profile")?;
    if r_e <= 0.0 || n <= 0.0 {
        return Err(BrahmandaError::InvalidGalaxy(
            "sersic_profile: r_e and n must be positive".to_string(),
        ));
    }
    // Approximation: b_n ≈ 2n - 1/3 + 4/(405n) (Ciotti & Bertin 1999)
    let b_n = 2.0 * n - 1.0 / 3.0 + 4.0 / (405.0 * n);
    let ratio = (r / r_e).powf(1.0 / n);
    ensure_finite((-b_n * (ratio - 1.0)).exp(), "sersic_profile")
}

/// Faber-Jackson relation: L ∝ σ⁴ (for elliptical galaxies).
///
/// Returns luminosity ratio L/L_ref = (σ/σ_ref)⁴.
///
/// # Arguments
/// * `sigma` — Velocity dispersion (km/s).
/// * `sigma_ref` — Reference velocity dispersion (km/s).
///
/// ```
/// use brahmanda::morphology::faber_jackson_ratio;
///
/// // Double the velocity dispersion → 16× the luminosity
/// let ratio = faber_jackson_ratio(400.0, 200.0).unwrap();
/// assert!((ratio - 16.0).abs() < 1e-6);
/// ```
#[inline]
pub fn faber_jackson_ratio(sigma: f64, sigma_ref: f64) -> Result<f64, BrahmandaError> {
    require_finite(sigma, "faber_jackson_ratio")?;
    require_finite(sigma_ref, "faber_jackson_ratio")?;
    if sigma <= 0.0 || sigma_ref <= 0.0 {
        return Err(BrahmandaError::InvalidGalaxy(
            "faber_jackson_ratio: dispersions must be positive".to_string(),
        ));
    }
    ensure_finite((sigma / sigma_ref).powi(4), "faber_jackson_ratio")
}

/// Tully-Fisher relation: L ∝ v_rot^α (for spiral galaxies).
///
/// Returns luminosity ratio L/L_ref = (v/v_ref)^alpha.
/// Typical α ≈ 4 (B-band) or α ≈ 3.5 (infrared).
///
/// ```
/// use brahmanda::morphology::tully_fisher_ratio;
///
/// let ratio = tully_fisher_ratio(300.0, 150.0, 4.0).unwrap();
/// assert!((ratio - 16.0).abs() < 1e-6);
/// ```
#[inline]
pub fn tully_fisher_ratio(
    v_rot: f64,
    v_ref: f64,
    alpha: f64,
) -> Result<f64, BrahmandaError> {
    require_finite(v_rot, "tully_fisher_ratio")?;
    require_finite(v_ref, "tully_fisher_ratio")?;
    require_finite(alpha, "tully_fisher_ratio")?;
    if v_rot <= 0.0 || v_ref <= 0.0 {
        return Err(BrahmandaError::InvalidGalaxy(
            "tully_fisher_ratio: velocities must be positive".to_string(),
        ));
    }
    ensure_finite((v_rot / v_ref).powf(alpha), "tully_fisher_ratio")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sersic_at_effective_radius() {
        // At R = R_e, profile should return ~1.0 (I/I_e = 1)
        let i = sersic_profile(1.0, 1.0, 4.0).unwrap();
        assert!((i - 1.0).abs() < 0.01, "Sersic at R_e: {i}");
    }

    #[test]
    fn test_sersic_decreases_outward() {
        let i1 = sersic_profile(1.0, 1.0, 4.0).unwrap();
        let i2 = sersic_profile(3.0, 1.0, 4.0).unwrap();
        assert!(i2 < i1);
    }

    #[test]
    fn test_faber_jackson_scaling() {
        // Double σ → 16× luminosity
        let ratio = faber_jackson_ratio(400.0, 200.0).unwrap();
        assert!((ratio - 16.0).abs() < 1e-6);
    }

    #[test]
    fn test_tully_fisher_scaling() {
        let ratio = tully_fisher_ratio(300.0, 150.0, 4.0).unwrap();
        assert!((ratio - 16.0).abs() < 1e-6);
    }

    #[test]
    fn test_invalid_inputs_rejected() {
        assert!(sersic_profile(1.0, -1.0, 4.0).is_err());
        assert!(sersic_profile(1.0, 1.0, 0.0).is_err());
        assert!(faber_jackson_ratio(-1.0, 200.0).is_err());
        assert!(tully_fisher_ratio(f64::NAN, 150.0, 4.0).is_err());
    }
}
