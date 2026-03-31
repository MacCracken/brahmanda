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

/// HSW void density profile — Hamaus, Sutter & Wandelt (2014).
///
/// Empirical model for the radial density contrast of cosmic voids:
///
/// δ(r) = δ_c × (1 - (r/r_s)^α) / (1 + (r/r_v)^β)
///
/// where r_s is the scale radius and r_v is the void radius.
/// The profile is underdense (δ < 0) in the interior and has a
/// compensating overdense ridge near the void edge.
///
/// # Arguments
/// * `r_mpc` — Distance from void center (Mpc).
/// * `r_v_mpc` — Void effective radius (Mpc).
/// * `delta_c` — Central density contrast (typically -0.8 to -0.4).
/// * `alpha` — Inner slope (typically ~2.0).
/// * `beta` — Outer steepness (typically ~6.0–10.0).
/// * `r_s_mpc` — Scale radius, transition point (typically ~0.9 × r_v).
///
/// ```
/// use brahmanda::cosmic_web::void_density_profile_hsw;
///
/// // Interior of a void: underdense
/// let d_center = void_density_profile_hsw(0.0, 20.0, -0.8, 2.0, 8.0, 18.0).unwrap();
/// assert!(d_center < 0.0);
///
/// // Far outside the void: approaches mean density (δ → 0)
/// let d_far = void_density_profile_hsw(100.0, 20.0, -0.8, 2.0, 8.0, 18.0).unwrap();
/// assert!(d_far.abs() < 0.1);
/// ```
pub fn void_density_profile_hsw(
    r_mpc: f64,
    r_v_mpc: f64,
    delta_c: f64,
    alpha: f64,
    beta: f64,
    r_s_mpc: f64,
) -> Result<f64, BrahmandaError> {
    require_finite(r_mpc, "void_density_profile_hsw")?;
    require_finite(r_v_mpc, "void_density_profile_hsw")?;
    require_finite(delta_c, "void_density_profile_hsw")?;
    require_finite(alpha, "void_density_profile_hsw")?;
    require_finite(beta, "void_density_profile_hsw")?;
    require_finite(r_s_mpc, "void_density_profile_hsw")?;
    if r_mpc < 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "void_density_profile_hsw: radius must be non-negative".to_string(),
        ));
    }
    if r_v_mpc <= 0.0 || r_s_mpc <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "void_density_profile_hsw: r_v and r_s must be positive".to_string(),
        ));
    }
    if alpha <= 0.0 || beta <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "void_density_profile_hsw: alpha and beta must be positive".to_string(),
        ));
    }

    let x_s = r_mpc / r_s_mpc;
    let x_v = r_mpc / r_v_mpc;
    let delta = delta_c * (1.0 - x_s.powf(alpha)) / (1.0 + x_v.powf(beta));
    ensure_finite(delta, "void_density_profile_hsw")
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

    #[test]
    fn test_hsw_center_underdense() {
        let d = void_density_profile_hsw(0.0, 20.0, -0.8, 2.0, 8.0, 18.0).unwrap();
        assert!((d - (-0.8)).abs() < 1e-10, "center should be δ_c: {d}");
    }

    #[test]
    fn test_hsw_far_field_approaches_zero() {
        let d = void_density_profile_hsw(200.0, 20.0, -0.8, 2.0, 8.0, 18.0).unwrap();
        assert!(d.abs() < 0.01, "far field should be near 0: {d}");
    }

    #[test]
    fn test_hsw_compensation_ridge() {
        // Near the void edge (r ~ r_s), the profile should cross zero
        // and become slightly positive (compensating overdensity)
        let d_inner = void_density_profile_hsw(10.0, 20.0, -0.8, 2.0, 8.0, 18.0).unwrap();
        let d_edge = void_density_profile_hsw(20.0, 20.0, -0.8, 2.0, 8.0, 18.0).unwrap();
        assert!(d_inner < d_edge, "density should increase toward edge");
    }

    #[test]
    fn test_hsw_invalid_inputs() {
        assert!(void_density_profile_hsw(-1.0, 20.0, -0.8, 2.0, 8.0, 18.0).is_err());
        assert!(void_density_profile_hsw(5.0, 0.0, -0.8, 2.0, 8.0, 18.0).is_err());
        assert!(void_density_profile_hsw(5.0, 20.0, -0.8, 0.0, 8.0, 18.0).is_err());
        assert!(void_density_profile_hsw(f64::NAN, 20.0, -0.8, 2.0, 8.0, 18.0).is_err());
    }
}
