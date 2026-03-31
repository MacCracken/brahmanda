//! Cosmic web topology — filaments, voids, sheets, nodes.
//!
//! Classification and statistical characterization of large-scale structure.

use serde::{Deserialize, Serialize};

use crate::error::{BrahmandaError, ensure_finite, require_all_finite, require_finite};

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

impl core::fmt::Display for WebEnvironment {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Node => f.write_str("Node"),
            Self::Filament => f.write_str("Filament"),
            Self::Sheet => f.write_str("Sheet"),
            Self::Void => f.write_str("Void"),
        }
    }
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
    ensure_finite((r0_mpc / r_mpc).powf(gamma), "two_point_correlation")
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

// ============================================================
// Void-in-void excursion set (Sheth & van de Weygaert 2004)
// ============================================================

/// Shell-crossing threshold for void formation (linear theory).
pub const DELTA_V: f64 = -2.717;

/// Void size distribution — Sheth & van de Weygaert (2004).
///
/// The void abundance in the excursion set framework with two barriers:
/// void shell-crossing at δ_v ≈ -2.717 and halo collapse at δ_c = 1.686.
///
/// Returns the volume fraction occupied by voids of radius R at redshift z,
/// expressed as dn/dlnR (comoving number density per unit ln R).
///
/// Uses the first-crossing distribution for the void barrier,
/// suppressed by the void-in-cloud correction:
///
/// f(S) = |δ_v|/√(2πS) × exp(-δ_v²/(2S)) × [1 - exp(-|δ_v|/δ_c × D²)]
///
/// where S = σ²(R, z) and D = (δ_c² - δ_v²)/(2S).
///
/// # Arguments
/// * `r_mpc` — Void radius in Mpc/h.
/// * `z` — Redshift.
///
/// ```
/// use brahmanda::cosmic_web::void_abundance_svdw;
///
/// let n = void_abundance_svdw(10.0, 0.0).unwrap();
/// assert!(n > 0.0, "void abundance should be positive");
///
/// // Larger voids are rarer
/// let n_big = void_abundance_svdw(30.0, 0.0).unwrap();
/// assert!(n_big < n);
/// ```
pub fn void_abundance_svdw(r_mpc: f64, z: f64) -> Result<f64, BrahmandaError> {
    require_finite(r_mpc, "void_abundance_svdw")?;
    require_finite(z, "void_abundance_svdw")?;
    if r_mpc <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "void_abundance_svdw: radius must be positive".to_string(),
        ));
    }

    let sigma = crate::power_spectrum::sigma_r(r_mpc, z)?;
    let s = sigma * sigma;
    let delta_c = crate::halo::DELTA_C;
    let dv = DELTA_V.abs();

    // First-crossing for void barrier (Gaussian)
    let f_void = dv / (2.0 * std::f64::consts::PI * s).sqrt() * (-dv * dv / (2.0 * s)).exp();

    // Void-in-cloud suppression: fraction of voids NOT embedded in a
    // collapsing region. The suppression factor is approximately
    // exp(-δ_c² / (2S)), the probability the trajectory stays below δ_c.
    let vic_suppression = (-delta_c * delta_c / (2.0 * s)).exp();

    ensure_finite(f_void * vic_suppression, "void_abundance_svdw")
}

// ============================================================
// Filament spine extraction (Hessian-based)
// ============================================================

/// Filament classification from density Hessian eigenvalues.
///
/// Classifies the local geometry based on the eigenvalues of the
/// density Hessian matrix (second derivatives of the density field).
///
/// Filaments are detected where one eigenvalue is much larger than
/// the other two (ridge condition).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum HessianMorphology {
    /// Blob/cluster: all eigenvalues large and negative (local max).
    Blob,
    /// Filament: one eigenvalue ~ 0, two large negative (ridge).
    Filament,
    /// Sheet/wall: two eigenvalues ~ 0, one large negative (plane).
    Wall,
    /// Void: all eigenvalues ~ 0 or positive (local min/saddle).
    Background,
}

impl core::fmt::Display for HessianMorphology {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Blob => f.write_str("Blob"),
            Self::Filament => f.write_str("Filament"),
            Self::Wall => f.write_str("Wall"),
            Self::Background => f.write_str("Background"),
        }
    }
}

/// Classify local morphology from Hessian eigenvalues.
///
/// Uses the ratio of eigenvalues to distinguish filaments from blobs,
/// walls, and background. Eigenvalues should be sorted: |λ₁| ≥ |λ₂| ≥ |λ₃|.
///
/// # Arguments
/// * `eigenvalues` — Three eigenvalues of the density Hessian, sorted by absolute value descending.
/// * `ridge_threshold` — Ratio threshold for ridge detection (typical: 0.5).
///
/// ```
/// use brahmanda::cosmic_web::{classify_hessian_morphology, HessianMorphology};
///
/// // Filament: one small eigenvalue, two large negative
/// let m = classify_hessian_morphology(&[-10.0, -8.0, -0.1], 0.5).unwrap();
/// assert_eq!(m, HessianMorphology::Filament);
///
/// // Blob: all eigenvalues large and negative
/// let m = classify_hessian_morphology(&[-10.0, -9.0, -8.0], 0.5).unwrap();
/// assert_eq!(m, HessianMorphology::Blob);
/// ```
pub fn classify_hessian_morphology(
    eigenvalues: &[f64; 3],
    ridge_threshold: f64,
) -> Result<HessianMorphology, BrahmandaError> {
    require_all_finite(eigenvalues, "classify_hessian_morphology")?;
    require_finite(ridge_threshold, "classify_hessian_morphology")?;
    if ridge_threshold <= 0.0 || ridge_threshold >= 1.0 {
        return Err(BrahmandaError::InvalidStructure(
            "classify_hessian_morphology: ridge_threshold must be in (0, 1)".to_string(),
        ));
    }

    let mut abs_ev = [
        eigenvalues[0].abs(),
        eigenvalues[1].abs(),
        eigenvalues[2].abs(),
    ];
    abs_ev.sort_by(|a, b| b.total_cmp(a));

    let lam_max = abs_ev[0];
    if lam_max < 1e-30 {
        return Ok(HessianMorphology::Background);
    }

    let r1 = abs_ev[1] / lam_max; // ratio of 2nd to 1st
    let r2 = abs_ev[2] / lam_max; // ratio of 3rd to 1st

    // All three eigenvalues are negative and comparable → blob
    let all_negative = eigenvalues.iter().all(|&e| e < 0.0);

    if all_negative && r1 > ridge_threshold && r2 > ridge_threshold {
        Ok(HessianMorphology::Blob)
    } else if all_negative && r1 > ridge_threshold && r2 <= ridge_threshold {
        // Two large negative, one small → filament
        Ok(HessianMorphology::Filament)
    } else if all_negative && r1 <= ridge_threshold {
        // One large negative, two small → wall
        Ok(HessianMorphology::Wall)
    } else {
        Ok(HessianMorphology::Background)
    }
}

/// Filamentarity parameter: F = (λ₂ - λ₃) / (λ₁ - λ₃).
///
/// Measures how filament-like the local geometry is.
/// F → 1 for perfect filaments, F → 0 for sheets or blobs.
///
/// # Arguments
/// * `eigenvalues` — Three Hessian eigenvalues sorted: λ₁ ≤ λ₂ ≤ λ₃ (most negative first).
///
/// ```
/// use brahmanda::cosmic_web::filamentarity;
///
/// // Strong filament: λ₁ ≈ λ₂ << λ₃
/// let f = filamentarity(&[-10.0, -9.0, -0.1]).unwrap();
/// assert!(f > 0.8);
///
/// // Isotropic (blob): λ₁ ≈ λ₂ ≈ λ₃
/// let f = filamentarity(&[-10.0, -10.0, -10.0]).unwrap();
/// assert!(f < 0.01);
/// ```
pub fn filamentarity(eigenvalues: &[f64; 3]) -> Result<f64, BrahmandaError> {
    require_all_finite(eigenvalues, "filamentarity")?;

    // Sort ascending (most negative first)
    let mut ev = *eigenvalues;
    ev.sort_by(|a, b| a.total_cmp(b));

    let denom = ev[0] - ev[2]; // λ₁ - λ₃ (negative for overdensities)
    if denom.abs() < 1e-30 {
        return Ok(0.0); // isotropic
    }
    ensure_finite((ev[1] - ev[2]) / denom, "filamentarity")
}

// ============================================================
// Minkowski functionals
// ============================================================

/// Minkowski functionals for a Gaussian random field excursion set.
///
/// For a 3D Gaussian field with RMS σ₀ and spectral moments σ₁, σ₂,
/// the four Minkowski functionals of the excursion set above threshold
/// ν = δ/σ₀ are analytically known (Tomita 1986, Matsubara 2003).
///
/// Returns (V₀, V₁, V₂, V₃):
/// - V₀ — Volume fraction
/// - V₁ — Surface area (per unit volume)
/// - V₂ — Integrated mean curvature (per unit volume)
/// - V₃ — Euler characteristic (per unit volume)
///
/// # Arguments
/// * `nu` — Threshold in units of σ₀ (ν = δ/σ₀).
/// * `sigma0` — RMS fluctuation σ₀.
/// * `sigma1` — First spectral moment σ₁ = ⟨k²⟩^(1/2) σ₀ (related to gradient).
///
/// ```
/// use brahmanda::cosmic_web::minkowski_functionals;
///
/// let (v0, v1, v2, v3) = minkowski_functionals(0.0, 1.0, 0.5).unwrap();
/// // At ν=0, half the volume is above threshold
/// assert!((v0 - 0.5).abs() < 0.01);
/// assert!(v1 > 0.0); // non-zero surface area
/// ```
pub fn minkowski_functionals(
    nu: f64,
    sigma0: f64,
    sigma1: f64,
) -> Result<(f64, f64, f64, f64), BrahmandaError> {
    require_finite(nu, "minkowski_functionals")?;
    require_finite(sigma0, "minkowski_functionals")?;
    require_finite(sigma1, "minkowski_functionals")?;
    if sigma0 <= 0.0 || sigma1 <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "minkowski_functionals: sigma0 and sigma1 must be positive".to_string(),
        ));
    }

    let pi = std::f64::consts::PI;
    let exp_term = (-nu * nu / 2.0).exp();
    let sqrt_2pi = (2.0 * pi).sqrt();

    // Spectral parameter λ = σ₁/σ₀
    let lambda = sigma1 / sigma0;

    // V₀: volume fraction = erfc(ν/√2) / 2
    // Use 1 - erf approximation: erfc(x) ≈ (2/√π) exp(-x²) / (x + √(x²+4/π)) for x > 0
    // For general ν, compute via the Gaussian CDF
    let v0 = 0.5 * erfc_approx(nu / std::f64::consts::SQRT_2);

    // V₁: surface area density
    let v1 = lambda / (6.0 * pi) * exp_term;

    // V₂: integrated mean curvature
    let v2 = lambda * lambda / (6.0 * pi * sqrt_2pi) * nu * exp_term;

    // V₃: Euler characteristic density
    let v3 = lambda.powi(3) / (6.0 * pi * pi) * (nu * nu - 1.0) * exp_term;

    Ok((
        ensure_finite(v0, "minkowski_functionals:V0")?,
        ensure_finite(v1, "minkowski_functionals:V1")?,
        ensure_finite(v2, "minkowski_functionals:V2")?,
        ensure_finite(v3, "minkowski_functionals:V3")?,
    ))
}

/// Complementary error function approximation.
///
/// erfc(x) = 1 - erf(x), using Abramowitz & Stegun 7.1.26.
fn erfc_approx(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc_approx(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    poly * (-x * x).exp()
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

    // -- void excursion set --

    #[test]
    fn test_void_abundance_positive() {
        let n = void_abundance_svdw(10.0, 0.0).unwrap();
        assert!(n > 0.0, "void abundance should be positive: {n}");
    }

    #[test]
    fn test_void_abundance_decreases_with_size() {
        let n1 = void_abundance_svdw(5.0, 0.0).unwrap();
        let n2 = void_abundance_svdw(15.0, 0.0).unwrap();
        let n3 = void_abundance_svdw(30.0, 0.0).unwrap();
        assert!(n1 > n2 && n2 > n3, "larger voids should be rarer");
    }

    #[test]
    fn test_void_abundance_invalid() {
        assert!(void_abundance_svdw(0.0, 0.0).is_err());
        assert!(void_abundance_svdw(-1.0, 0.0).is_err());
        assert!(void_abundance_svdw(f64::NAN, 0.0).is_err());
    }

    // -- hessian morphology --

    #[test]
    fn test_hessian_blob() {
        let m = classify_hessian_morphology(&[-10.0, -9.0, -8.0], 0.5).unwrap();
        assert_eq!(m, HessianMorphology::Blob);
    }

    #[test]
    fn test_hessian_filament() {
        let m = classify_hessian_morphology(&[-10.0, -8.0, -0.1], 0.5).unwrap();
        assert_eq!(m, HessianMorphology::Filament);
    }

    #[test]
    fn test_hessian_wall() {
        let m = classify_hessian_morphology(&[-10.0, -0.1, -0.1], 0.5).unwrap();
        assert_eq!(m, HessianMorphology::Wall);
    }

    #[test]
    fn test_hessian_background() {
        let m = classify_hessian_morphology(&[1.0, 0.5, 0.1], 0.5).unwrap();
        assert_eq!(m, HessianMorphology::Background);
    }

    #[test]
    fn test_hessian_invalid() {
        assert!(classify_hessian_morphology(&[f64::NAN, 0.0, 0.0], 0.5).is_err());
        assert!(classify_hessian_morphology(&[1.0, 0.0, 0.0], 0.0).is_err());
        assert!(classify_hessian_morphology(&[1.0, 0.0, 0.0], 1.0).is_err());
    }

    #[test]
    fn test_filamentarity_strong() {
        let f = filamentarity(&[-10.0, -9.0, -0.1]).unwrap();
        assert!(f > 0.8, "strong filament: {f}");
    }

    #[test]
    fn test_filamentarity_isotropic() {
        let f = filamentarity(&[-5.0, -5.0, -5.0]).unwrap();
        assert!(f.abs() < 0.01, "isotropic: {f}");
    }

    #[test]
    fn test_filamentarity_sheet() {
        // Sheet: λ₁ << λ₂ ≈ λ₃
        let f = filamentarity(&[-10.0, -0.1, -0.1]).unwrap();
        assert!(f < 0.02, "sheet should have low filamentarity: {f}");
    }

    // -- minkowski functionals --

    #[test]
    fn test_minkowski_v0_at_zero() {
        let (v0, _, _, _) = minkowski_functionals(0.0, 1.0, 0.5).unwrap();
        assert!((v0 - 0.5).abs() < 0.01, "V₀(ν=0) ≈ 0.5: {v0}");
    }

    #[test]
    fn test_minkowski_v0_monotonic() {
        let (v0_neg, _, _, _) = minkowski_functionals(-1.0, 1.0, 0.5).unwrap();
        let (v0_zero, _, _, _) = minkowski_functionals(0.0, 1.0, 0.5).unwrap();
        let (v0_pos, _, _, _) = minkowski_functionals(1.0, 1.0, 0.5).unwrap();
        assert!(
            v0_neg > v0_zero && v0_zero > v0_pos,
            "V₀ should decrease with ν"
        );
    }

    #[test]
    fn test_minkowski_v3_genus() {
        // V₃ ∝ (ν²-1), so V₃(ν=0) < 0, V₃(ν=±√2) ≈ 0
        let (_, _, _, v3_0) = minkowski_functionals(0.0, 1.0, 0.5).unwrap();
        assert!(v3_0 < 0.0, "V₃(ν=0) should be negative: {v3_0}");

        let (_, _, _, v3_high) = minkowski_functionals(2.0, 1.0, 0.5).unwrap();
        assert!(v3_high > 0.0, "V₃(ν=2) should be positive: {v3_high}");
    }

    #[test]
    fn test_minkowski_invalid() {
        assert!(minkowski_functionals(0.0, 0.0, 0.5).is_err());
        assert!(minkowski_functionals(0.0, 1.0, -1.0).is_err());
        assert!(minkowski_functionals(f64::NAN, 1.0, 0.5).is_err());
    }
}
