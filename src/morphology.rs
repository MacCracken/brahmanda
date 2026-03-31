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

/// Mass-metallicity relation — Tremonti et al. (2004).
///
/// 12 + log₁₀(O/H) as a function of stellar mass, with an asymptotic
/// turnover at high mass. Calibrated for z ≈ 0 star-forming galaxies.
///
/// Uses the polynomial fit:
/// 12 + log(O/H) = -1.492 + 1.847 x - 0.08026 x²
/// where x = log₁₀(M* / M_sun).
///
/// # Arguments
/// * `log_m_star` — log₁₀(M*/M_sun), valid range ~8–12.
///
/// ```
/// use brahmanda::morphology::mass_metallicity_tremonti;
///
/// // Solar neighborhood: M* ~ 10^10.5 → 12+log(O/H) ≈ 8.9
/// let z_met = mass_metallicity_tremonti(10.5).unwrap();
/// assert!(z_met > 8.5 && z_met < 9.3, "MZR at 10.5: {z_met}");
///
/// // Dwarf galaxy: lower metallicity
/// let z_dwarf = mass_metallicity_tremonti(8.0).unwrap();
/// assert!(z_dwarf < z_met);
/// ```
pub fn mass_metallicity_tremonti(log_m_star: f64) -> Result<f64, BrahmandaError> {
    require_finite(log_m_star, "mass_metallicity_tremonti")?;
    let x = log_m_star;
    let oh = -1.492 + 1.847 * x - 0.08026 * x * x;
    ensure_finite(oh, "mass_metallicity_tremonti")
}

/// Mass-metallicity relation with redshift evolution.
///
/// Extends the Tremonti (2004) z=0 relation to higher redshift using
/// the observed ~0.3 dex offset per unit redshift (Zahid et al. 2014):
///
/// 12 + log(O/H)(z) ≈ MZR(z=0) - 0.3 × z
///
/// ```
/// use brahmanda::morphology::mass_metallicity_z;
///
/// let z0 = mass_metallicity_z(10.5, 0.0).unwrap();
/// let z1 = mass_metallicity_z(10.5, 1.0).unwrap();
/// assert!(z0 > z1, "higher z → lower metallicity");
/// assert!((z0 - z1 - 0.3).abs() < 0.01);
/// ```
pub fn mass_metallicity_z(log_m_star: f64, z: f64) -> Result<f64, BrahmandaError> {
    require_finite(z, "mass_metallicity_z")?;
    if z < 0.0 {
        return Err(BrahmandaError::InvalidGalaxy(
            "mass_metallicity_z: z must be non-negative".to_string(),
        ));
    }
    let mzr_0 = mass_metallicity_tremonti(log_m_star)?;
    ensure_finite(mzr_0 - 0.3 * z, "mass_metallicity_z")
}

/// Cosmic star formation rate density — Madau & Dickinson (2014).
///
/// ρ_SFR(z) = 0.015 × (1+z)^2.7 / [1 + ((1+z)/2.9)^5.6]
///
/// Returns the star formation rate density in M_sun yr⁻¹ Mpc⁻³.
///
/// ```
/// use brahmanda::morphology::sfr_density_madau_dickinson;
///
/// // SFR peaks at z ≈ 2
/// let sfr_0 = sfr_density_madau_dickinson(0.0).unwrap();
/// let sfr_2 = sfr_density_madau_dickinson(2.0).unwrap();
/// let sfr_5 = sfr_density_madau_dickinson(5.0).unwrap();
/// assert!(sfr_2 > sfr_0, "SFR peaks around z~2");
/// assert!(sfr_2 > sfr_5, "SFR declines at high z");
/// ```
pub fn sfr_density_madau_dickinson(z: f64) -> Result<f64, BrahmandaError> {
    require_finite(z, "sfr_density_madau_dickinson")?;
    if z < 0.0 {
        return Err(BrahmandaError::InvalidGalaxy(
            "sfr_density_madau_dickinson: z must be non-negative".to_string(),
        ));
    }
    let opz = 1.0 + z;
    let rho_sfr = 0.015 * opz.powf(2.7) / (1.0 + (opz / 2.9).powf(5.6));
    ensure_finite(rho_sfr, "sfr_density_madau_dickinson")
}

/// Cumulative stellar mass density — integrated Madau-Dickinson SFR.
///
/// ρ*(z) = ∫_z^∞ ρ_SFR(z') |dt/dz'| dz'
///
/// Approximate: integrates from z to z_max=10 using the SFR density
/// and cosmic time relation |dt/dz| = 1/((1+z)H(z)).
///
/// Returns stellar mass density in M_sun Mpc⁻³.
///
/// ```
/// use brahmanda::morphology::stellar_mass_density;
///
/// let rho_0 = stellar_mass_density(0.0).unwrap();
/// let rho_2 = stellar_mass_density(2.0).unwrap();
/// assert!(rho_0 > rho_2, "more stars at z=0 than z=2");
/// assert!(rho_0 > 1e7, "ρ*(z=0) should be substantial");
/// ```
pub fn stellar_mass_density(z: f64) -> Result<f64, BrahmandaError> {
    require_finite(z, "stellar_mass_density")?;
    if z < 0.0 {
        return Err(BrahmandaError::InvalidGalaxy(
            "stellar_mass_density: z must be non-negative".to_string(),
        ));
    }

    let z_max = 10.0;
    if z >= z_max {
        return Ok(0.0);
    }

    let n = 200_usize;
    let dz = (z_max - z) / n as f64;
    let mut sum = 0.0;

    // |dt/dz| = 1/((1+z)H(z)), H(z) = H₀ E(z)
    // H₀ in yr⁻¹: H₀ = 67.4 km/s/Mpc = 67.4/(3.0857e19 km) s⁻¹ → × 3.156e7 s/yr⁻¹
    let h0_per_yr = crate::constants::H0_SI * 3.15576e7; // s⁻¹ → yr⁻¹

    for i in 0..n {
        let z0 = z + i as f64 * dz;
        let z1 = z0 + dz / 2.0;
        let z2 = z0 + dz;

        let integrand = |zz: f64| -> Result<f64, BrahmandaError> {
            let sfr = sfr_density_madau_dickinson(zz)?;
            let e = crate::power_spectrum::hubble_parameter_ratio(
                zz, crate::constants::OMEGA_M, -1.0, 0.0,
            )?;
            Ok(sfr / ((1.0 + zz) * h0_per_yr * e))
        };

        let f0 = integrand(z0)?;
        let f1 = integrand(z1)?;
        let f2 = integrand(z2)?;
        sum += (dz / 6.0) * (f0 + 4.0 * f1 + f2);
    }

    ensure_finite(sum, "stellar_mass_density")
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

    // -- mass-metallicity --

    #[test]
    fn test_mzr_increases_with_mass() {
        let z8 = mass_metallicity_tremonti(8.0).unwrap();
        let z10 = mass_metallicity_tremonti(10.0).unwrap();
        let z11 = mass_metallicity_tremonti(11.0).unwrap();
        // MZR has turnover, but should increase from 8 → 10
        assert!(z10 > z8, "metallicity should increase with mass");
        // Near turnover at high mass
        assert!(z11 > z8);
    }

    #[test]
    fn test_mzr_solar_neighborhood() {
        let z = mass_metallicity_tremonti(10.5).unwrap();
        assert!(z > 8.5 && z < 9.3, "solar neighborhood: {z}");
    }

    #[test]
    fn test_mzr_z_evolution() {
        let z0 = mass_metallicity_z(10.0, 0.0).unwrap();
        let z1 = mass_metallicity_z(10.0, 1.0).unwrap();
        let z2 = mass_metallicity_z(10.0, 2.0).unwrap();
        assert!(z0 > z1 && z1 > z2, "metallicity decreases with z");
    }

    #[test]
    fn test_mzr_invalid() {
        assert!(mass_metallicity_tremonti(f64::NAN).is_err());
        assert!(mass_metallicity_z(10.0, -1.0).is_err());
    }

    // -- SFR density --

    #[test]
    fn test_sfr_peaks_around_z2() {
        let sfr_0 = sfr_density_madau_dickinson(0.0).unwrap();
        let sfr_2 = sfr_density_madau_dickinson(2.0).unwrap();
        let sfr_5 = sfr_density_madau_dickinson(5.0).unwrap();
        assert!(sfr_2 > sfr_0 && sfr_2 > sfr_5, "SFR should peak near z~2");
    }

    #[test]
    fn test_sfr_z0_value() {
        let sfr = sfr_density_madau_dickinson(0.0).unwrap();
        // ρ_SFR(z=0) ≈ 0.015 M_sun/yr/Mpc³
        assert!(sfr > 0.005 && sfr < 0.03, "ρ_SFR(z=0) = {sfr}");
    }

    #[test]
    fn test_sfr_positive() {
        for z in [0.0, 1.0, 3.0, 6.0, 10.0] {
            let sfr = sfr_density_madau_dickinson(z).unwrap();
            assert!(sfr > 0.0);
        }
    }

    #[test]
    fn test_sfr_invalid() {
        assert!(sfr_density_madau_dickinson(-1.0).is_err());
        assert!(sfr_density_madau_dickinson(f64::NAN).is_err());
    }

    // -- stellar mass density --

    #[test]
    fn test_stellar_mass_density_decreases_with_z() {
        let rho_0 = stellar_mass_density(0.0).unwrap();
        let rho_1 = stellar_mass_density(1.0).unwrap();
        let rho_3 = stellar_mass_density(3.0).unwrap();
        assert!(rho_0 > rho_1 && rho_1 > rho_3);
    }

    #[test]
    fn test_stellar_mass_density_positive() {
        let rho = stellar_mass_density(0.0).unwrap();
        assert!(rho > 0.0);
    }

    #[test]
    fn test_stellar_mass_density_invalid() {
        assert!(stellar_mass_density(-1.0).is_err());
    }
}
