//! Cosmic web topology — filaments, voids, sheets, nodes.
//!
//! Classification and statistical characterization of large-scale structure.

use serde::{Deserialize, Serialize};

use crate::error::{ensure_finite, require_finite, BrahmandaError};

/// Classification of cosmic web environments.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum WebEnvironment {
    /// Dense intersection of filaments; galaxy clusters.
    Node,
    /// Elongated overdensity connecting nodes.
    Filament,
    /// Two-dimensional overdense surface.
    Sheet,
    /// Underdense region between filaments.
    Void,
}

/// Classify environment from eigenvalues of the tidal tensor.
///
/// The tidal (deformation) tensor has 3 eigenvalues (λ₁ ≥ λ₂ ≥ λ₃).
/// Classification by counting eigenvalues above threshold δ_th:
/// - 3 above → Node
/// - 2 above → Filament
/// - 1 above → Sheet
/// - 0 above → Void
///
/// # Arguments
/// * `lambda` — Three eigenvalues of the tidal tensor [λ₁, λ₂, λ₃] (order does not matter).
/// * `threshold` — Classification threshold δ_th (typically 0.0 or density-dependent).
///
/// ```
/// use brahmanda::cosmic_web::{classify_web_environment, WebEnvironment};
///
/// let env = classify_web_environment(&[2.0, 1.0, 0.5], 0.0).unwrap();
/// assert_eq!(env, WebEnvironment::Node);
///
/// let env = classify_web_environment(&[-0.1, -0.5, -0.8], 0.0).unwrap();
/// assert_eq!(env, WebEnvironment::Void);
/// ```
pub fn classify_web_environment(
    lambda: &[f64; 3],
    threshold: f64,
) -> Result<WebEnvironment, BrahmandaError> {
    for &l in lambda {
        require_finite(l, "classify_web_environment")?;
    }
    require_finite(threshold, "classify_web_environment")?;

    let count = lambda.iter().filter(|&&l| l > threshold).count();
    Ok(match count {
        3 => WebEnvironment::Node,
        2 => WebEnvironment::Filament,
        1 => WebEnvironment::Sheet,
        _ => WebEnvironment::Void,
    })
}

/// Void radius from void volume (Mpc).
///
/// Assumes spherical void: R = (3V / 4π)^(1/3).
///
/// ```
/// use brahmanda::cosmic_web::void_radius_mpc;
///
/// let r = void_radius_mpc(1000.0).unwrap();
/// assert!(r > 5.0 && r < 7.0); // ~6.2 Mpc
/// ```
#[inline]
pub fn void_radius_mpc(volume_mpc3: f64) -> Result<f64, BrahmandaError> {
    require_finite(volume_mpc3, "void_radius_mpc")?;
    if volume_mpc3 <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "void_radius_mpc: volume must be positive".to_string(),
        ));
    }
    ensure_finite(
        (3.0 * volume_mpc3 / (4.0 * std::f64::consts::PI)).cbrt(),
        "void_radius_mpc",
    )
}

/// Density contrast δ = (ρ - ρ̄) / ρ̄.
///
/// δ > 0 for overdense regions (nodes, filaments).
/// δ < 0 for underdense regions (voids).
/// δ = 0 for mean density.
///
/// ```
/// use brahmanda::cosmic_web::density_contrast;
///
/// let delta = density_contrast(2.0, 1.0).unwrap();
/// assert!((delta - 1.0).abs() < 1e-10); // overdense
///
/// let delta = density_contrast(0.5, 1.0).unwrap();
/// assert!((delta - (-0.5)).abs() < 1e-10); // underdense
/// ```
#[inline]
pub fn density_contrast(rho: f64, rho_mean: f64) -> Result<f64, BrahmandaError> {
    require_finite(rho, "density_contrast")?;
    require_finite(rho_mean, "density_contrast")?;
    if rho_mean <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "density_contrast: mean density must be positive".to_string(),
        ));
    }
    ensure_finite((rho - rho_mean) / rho_mean, "density_contrast")
}

/// Two-point correlation function ξ(r) — power-law approximation.
///
/// ξ(r) = (r₀/r)^γ where r₀ ≈ 5 h⁻¹ Mpc, γ ≈ 1.8 for galaxies.
///
/// Returns the excess probability of finding a galaxy pair at separation r
/// relative to a uniform distribution.
///
/// ```
/// use brahmanda::cosmic_web::two_point_correlation_power_law;
///
/// // At r = r₀, ξ = 1
/// let xi = two_point_correlation_power_law(5.0, 5.0, 1.8).unwrap();
/// assert!((xi - 1.0).abs() < 1e-10);
///
/// // Correlation decreases with distance
/// let xi_far = two_point_correlation_power_law(10.0, 5.0, 1.8).unwrap();
/// assert!(xi_far < xi);
/// ```
#[inline]
pub fn two_point_correlation_power_law(
    r_mpc: f64,
    r0_mpc: f64,
    gamma: f64,
) -> Result<f64, BrahmandaError> {
    require_finite(r_mpc, "two_point_correlation")?;
    require_finite(r0_mpc, "two_point_correlation")?;
    require_finite(gamma, "two_point_correlation")?;
    if r_mpc <= 0.0 || r0_mpc <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "two_point_correlation: distances must be positive".to_string(),
        ));
    }
    ensure_finite(
        (r0_mpc / r_mpc).powf(gamma),
        "two_point_correlation",
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_classify_node() {
        let env = classify_web_environment(&[1.0, 0.5, 0.3], 0.0).unwrap();
        assert_eq!(env, WebEnvironment::Node);
    }

    #[test]
    fn test_classify_filament() {
        let env = classify_web_environment(&[1.0, 0.5, -0.3], 0.0).unwrap();
        assert_eq!(env, WebEnvironment::Filament);
    }

    #[test]
    fn test_classify_sheet() {
        let env = classify_web_environment(&[1.0, -0.5, -0.3], 0.0).unwrap();
        assert_eq!(env, WebEnvironment::Sheet);
    }

    #[test]
    fn test_classify_void() {
        let env = classify_web_environment(&[-0.1, -0.5, -0.8], 0.0).unwrap();
        assert_eq!(env, WebEnvironment::Void);
    }

    #[test]
    fn test_density_contrast_overdense() {
        let delta = density_contrast(2.0, 1.0).unwrap();
        assert!((delta - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_density_contrast_void() {
        let delta = density_contrast(0.5, 1.0).unwrap();
        assert!((delta - (-0.5)).abs() < 1e-10);
    }

    #[test]
    fn test_correlation_decreases_with_distance() {
        let xi1 = two_point_correlation_power_law(5.0, 5.0, 1.8).unwrap();
        let xi2 = two_point_correlation_power_law(10.0, 5.0, 1.8).unwrap();
        assert!(xi1 > xi2);
    }

    #[test]
    fn test_correlation_at_r0() {
        // At r = r₀: ξ = 1
        let xi = two_point_correlation_power_law(5.0, 5.0, 1.8).unwrap();
        assert!((xi - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_void_radius() {
        // 1000 Mpc³ → R ≈ 6.2 Mpc
        let r = void_radius_mpc(1000.0).unwrap();
        assert!(r > 5.0 && r < 7.0);
    }
}
