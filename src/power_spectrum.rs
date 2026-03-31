//! Matter power spectrum — structure formation, transfer functions.
//!
//! The power spectrum P(k) describes the amplitude of density fluctuations
//! as a function of wavenumber k. Combined with the growth factor, it
//! predicts structure at all scales and redshifts.

use std::f64::consts::PI;

use crate::constants::{N_S, OMEGA_B, OMEGA_M, SIGMA_8};
use crate::error::{ensure_finite, require_finite, BrahmandaError};

/// Eisenstein-Hu transfer function T(k) — no-wiggle approximation.
///
/// Approximates the full transfer function including baryon suppression
/// but without BAO wiggles. Adequate for most structure formation work.
///
/// # Arguments
/// * `k_mpc` — Wavenumber in h/Mpc.
/// * `omega_m` — Matter density parameter.
/// * `omega_b` — Baryon density parameter.
/// * `h` — Dimensionless Hubble parameter (H₀/100).
pub fn transfer_function_eh(
    k_mpc: f64,
    omega_m: f64,
    omega_b: f64,
    h: f64,
) -> Result<f64, BrahmandaError> {
    require_finite(k_mpc, "transfer_function_eh")?;
    require_finite(omega_m, "transfer_function_eh")?;
    require_finite(omega_b, "transfer_function_eh")?;
    require_finite(h, "transfer_function_eh")?;
    if k_mpc <= 0.0 || omega_m <= 0.0 || h <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "transfer_function_eh: k, omega_m, h must be positive".to_string(),
        ));
    }

    let theta = 2.725 / 2.7; // CMB temperature ratio
    let omega_mh2 = omega_m * h * h;
    let omega_bh2 = omega_b * h * h;

    // Sound horizon (Eq. 26 of Eisenstein & Hu 1998)
    let s = 44.5 * (omega_mh2).sqrt().recip()
        * (1.0 + (omega_bh2 / 0.024).powf(0.46)).recip()
        * (omega_mh2).powf(0.25)
        * 100.0; // rough approximation in Mpc/h

    // Silk damping scale
    let alpha_gamma = 1.0 - 0.328 * (omega_bh2 / omega_mh2).ln() * (omega_b / omega_m)
        + 0.38 * (omega_bh2 / omega_mh2).ln() * (omega_b / omega_m).powi(2);

    // Shape parameter
    let gamma_eff = omega_m * h * (alpha_gamma + (1.0 - alpha_gamma) / (1.0 + (0.43 * k_mpc * s).powi(4)));

    // BBKS transfer function with effective shape
    let q = k_mpc * theta * theta / gamma_eff;
    let t = (1.0 + 3.89 * q + (16.1 * q).powi(2) + (5.46 * q).powi(3) + (6.71 * q).powi(4))
        .powf(-0.25)
        * (2.34 * q).ln_1p()
        / (2.34 * q);

    ensure_finite(t, "transfer_function_eh")
}

/// Primordial power spectrum: P_prim(k) ∝ k^n_s.
///
/// Returns the unnormalized primordial power at wavenumber k.
#[inline]
pub fn primordial_power(k_mpc: f64, n_s: f64) -> Result<f64, BrahmandaError> {
    require_finite(k_mpc, "primordial_power")?;
    require_finite(n_s, "primordial_power")?;
    if k_mpc <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "primordial_power: k must be positive".to_string(),
        ));
    }
    ensure_finite(k_mpc.powf(n_s), "primordial_power")
}

/// Linear growth factor D(z) — Carroll, Press & Turner (1992) approximation.
///
/// D(z) / D(0) for a flat ΛCDM cosmology.
///
/// # Arguments
/// * `z` — Redshift.
/// * `omega_m` — Matter density parameter (at z=0).
pub fn growth_factor(z: f64, omega_m: f64) -> Result<f64, BrahmandaError> {
    require_finite(z, "growth_factor")?;
    require_finite(omega_m, "growth_factor")?;
    if z < -1.0 {
        return Err(BrahmandaError::InvalidStructure(
            "growth_factor: z must be >= -1".to_string(),
        ));
    }

    let omega_lambda = 1.0 - omega_m;

    let growth_at = |zz: f64| -> f64 {
        let a = 1.0 / (1.0 + zz);
        let om_z = omega_m / (omega_m + omega_lambda * a * a * a);
        let ol_z = omega_lambda * a * a * a / (omega_m + omega_lambda * a * a * a);
        // CPT92 approximation
        let d = 2.5 * om_z / (om_z.powf(4.0 / 7.0) - ol_z + (1.0 + om_z / 2.0) * (1.0 + ol_z / 70.0));
        d * a
    };

    let d_z = growth_at(z);
    let d_0 = growth_at(0.0);
    ensure_finite(d_z / d_0, "growth_factor")
}

/// σ(R) — RMS mass fluctuation in a sphere of radius R (Mpc/h).
///
/// Simplified: scales σ₈ by the growth factor and a radius-dependent factor.
/// σ(R, z) ≈ σ₈ × D(z) × (8/R)^((n_s+3)/6)
///
/// This is an approximation; the exact computation requires integration
/// over the full power spectrum with a top-hat window function.
pub fn sigma_r(r_mpc: f64, z: f64) -> Result<f64, BrahmandaError> {
    require_finite(r_mpc, "sigma_r")?;
    if r_mpc <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "sigma_r: radius must be positive".to_string(),
        ));
    }
    let d = growth_factor(z, OMEGA_M)?;
    let scaling = (8.0 / r_mpc).powf((N_S + 3.0) / 6.0);
    ensure_finite(SIGMA_8 * d * scaling, "sigma_r")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_transfer_function_small_k() {
        // At very small k (large scales), T(k) → 1
        let t = transfer_function_eh(0.001, OMEGA_M, OMEGA_B, 0.674).unwrap();
        assert!(t > 0.8, "T(k→0) should approach 1, got {t}");
    }

    #[test]
    fn test_transfer_function_large_k() {
        // At large k (small scales), T(k) < 1 (suppressed)
        let t = transfer_function_eh(1.0, OMEGA_M, OMEGA_B, 0.674).unwrap();
        assert!(t < 0.5, "T(k=1) should be suppressed, got {t}");
    }

    #[test]
    fn test_growth_factor_z0() {
        // D(z=0) / D(0) = 1 by definition
        let d = growth_factor(0.0, OMEGA_M).unwrap();
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_growth_factor_decreases_with_z() {
        let d0 = growth_factor(0.0, OMEGA_M).unwrap();
        let d1 = growth_factor(1.0, OMEGA_M).unwrap();
        let d5 = growth_factor(5.0, OMEGA_M).unwrap();
        assert!(d0 > d1);
        assert!(d1 > d5);
    }

    #[test]
    fn test_sigma_8_at_z0() {
        // σ(R=8, z=0) should be σ₈
        let s = sigma_r(8.0, 0.0).unwrap();
        assert!((s - SIGMA_8).abs() < 0.01, "σ(8,0) = {s}, expected {SIGMA_8}");
    }

    #[test]
    fn test_sigma_decreases_with_z() {
        let s0 = sigma_r(8.0, 0.0).unwrap();
        let s1 = sigma_r(8.0, 1.0).unwrap();
        assert!(s0 > s1, "fluctuations should decrease at higher z");
    }

    #[test]
    fn test_invalid_rejected() {
        assert!(transfer_function_eh(-1.0, OMEGA_M, OMEGA_B, 0.674).is_err());
        assert!(growth_factor(-2.0, OMEGA_M).is_err());
        assert!(sigma_r(0.0, 0.0).is_err());
    }
}
