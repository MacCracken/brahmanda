//! Cosmological parameter struct and presets.
//!
//! The [`Cosmology`] struct encapsulates all cosmological parameters needed
//! by brahmanda's functions. Use presets like [`Cosmology::planck2018()`] or
//! construct a custom cosmology with [`Cosmology::new()`].

use std::f64::consts::PI;
use std::sync::OnceLock;

use serde::{Deserialize, Serialize};

use crate::constants::{G, MPC_KM};
use crate::error::BrahmandaError;

/// A complete set of cosmological parameters.
///
/// Contains both primary parameters (set by the user or a preset) and
/// derived quantities computed at construction time.
///
/// # Presets
///
/// ```
/// use brahmanda::cosmology::Cosmology;
///
/// let cosmo = Cosmology::planck2018();
/// assert!((cosmo.h - 0.674).abs() < 1e-10);
/// assert!((cosmo.omega_m - 0.315).abs() < 1e-10);
/// ```
///
/// # Custom cosmology
///
/// ```
/// use brahmanda::cosmology::Cosmology;
///
/// let cosmo = Cosmology::new(
///     0.7, 0.3, 0.045, 0.0, 9e-5, 0.96, 0.8, -1.0, 0.0, 2.725,
/// ).unwrap();
/// assert!(cosmo.omega_lambda > 0.69);
/// ```
#[derive(Serialize, Deserialize)]
pub struct Cosmology {
    /// Dimensionless Hubble parameter h = H_0 / (100 km/s/Mpc).
    pub h: f64,
    /// Total matter density parameter (CDM + baryons).
    pub omega_m: f64,
    /// Baryon density parameter.
    pub omega_b: f64,
    /// Curvature density parameter.
    pub omega_k: f64,
    /// Radiation density parameter.
    pub omega_r: f64,
    /// Scalar spectral index.
    pub n_s: f64,
    /// RMS mass fluctuations in 8 Mpc/h spheres.
    pub sigma_8: f64,
    /// Dark energy equation of state at z=0.
    pub w0: f64,
    /// Dark energy equation of state evolution parameter.
    pub wa: f64,
    /// CMB temperature today (K).
    pub t_cmb: f64,

    // -- derived --
    /// Hubble constant in SI (s^-1).
    pub h0_si: f64,
    /// Critical density today (kg/m^3).
    pub rho_crit: f64,
    /// Dark energy density parameter (derived from flatness).
    pub omega_lambda: f64,

    // -- cached --
    /// Power spectrum amplitude A_s, lazy-computed to match sigma_8.
    #[serde(skip)]
    a_s_cache: OnceLock<f64>,
}

impl core::fmt::Debug for Cosmology {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("Cosmology")
            .field("h", &self.h)
            .field("omega_m", &self.omega_m)
            .field("omega_b", &self.omega_b)
            .field("omega_k", &self.omega_k)
            .field("omega_r", &self.omega_r)
            .field("n_s", &self.n_s)
            .field("sigma_8", &self.sigma_8)
            .field("w0", &self.w0)
            .field("wa", &self.wa)
            .field("t_cmb", &self.t_cmb)
            .field("h0_si", &self.h0_si)
            .field("rho_crit", &self.rho_crit)
            .field("omega_lambda", &self.omega_lambda)
            .finish()
    }
}

impl Clone for Cosmology {
    fn clone(&self) -> Self {
        Self {
            h: self.h,
            omega_m: self.omega_m,
            omega_b: self.omega_b,
            omega_k: self.omega_k,
            omega_r: self.omega_r,
            n_s: self.n_s,
            sigma_8: self.sigma_8,
            w0: self.w0,
            wa: self.wa,
            t_cmb: self.t_cmb,
            h0_si: self.h0_si,
            rho_crit: self.rho_crit,
            omega_lambda: self.omega_lambda,
            a_s_cache: OnceLock::new(),
        }
    }
}

impl PartialEq for Cosmology {
    fn eq(&self, other: &Self) -> bool {
        self.h == other.h
            && self.omega_m == other.omega_m
            && self.omega_b == other.omega_b
            && self.omega_k == other.omega_k
            && self.omega_r == other.omega_r
            && self.n_s == other.n_s
            && self.sigma_8 == other.sigma_8
            && self.w0 == other.w0
            && self.wa == other.wa
            && self.t_cmb == other.t_cmb
    }
}

/// Window (filter) functions for smoothing the density field.
///
/// Used in the computation of σ(R) to select the smoothing kernel.
///
/// ```
/// use brahmanda::cosmology::FilterFunction;
///
/// let w = FilterFunction::TopHat.window(0.0);
/// assert!((w - 1.0).abs() < 1e-10);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum FilterFunction {
    /// Real-space top-hat: W(x) = 3(sin x - x cos x) / x^3.
    TopHat,
    /// Gaussian: W(x) = exp(-x^2/2).
    Gaussian,
    /// Sharp k-space filter: W(x) = 1 if x < 1, else 0.
    SharpK,
}

impl FilterFunction {
    /// Evaluate the window function at x = kR.
    ///
    /// ```
    /// use brahmanda::cosmology::FilterFunction;
    ///
    /// // Top-hat at x=0 → 1
    /// assert!((FilterFunction::TopHat.window(0.0) - 1.0).abs() < 1e-10);
    ///
    /// // Gaussian at x=0 → 1
    /// assert!((FilterFunction::Gaussian.window(0.0) - 1.0).abs() < 1e-10);
    /// ```
    #[inline]
    #[must_use]
    pub fn window(self, x: f64) -> f64 {
        match self {
            Self::TopHat => {
                if x.abs() < 1e-6 {
                    1.0 - x * x / 10.0
                } else {
                    3.0 * (x.sin() - x * x.cos()) / (x * x * x)
                }
            }
            Self::Gaussian => (-x * x / 2.0).exp(),
            Self::SharpK => {
                if x < 1.0 {
                    1.0
                } else {
                    0.0
                }
            }
        }
    }
}

impl Cosmology {
    /// Create a new cosmology with the given parameters.
    ///
    /// # Arguments
    /// * `h` — Dimensionless Hubble parameter (must be > 0).
    /// * `omega_m` — Matter density parameter (must be >= 0).
    /// * `omega_b` — Baryon density parameter (must be >= 0 and <= omega_m).
    /// * `omega_k` — Curvature density parameter.
    /// * `omega_r` — Radiation density parameter (must be >= 0).
    /// * `n_s` — Scalar spectral index.
    /// * `sigma_8` — RMS fluctuations at 8 Mpc/h (must be > 0).
    /// * `w0` — Dark energy EoS at z=0.
    /// * `wa` — Dark energy EoS evolution.
    /// * `t_cmb` — CMB temperature (K, must be > 0).
    ///
    /// ```
    /// use brahmanda::cosmology::Cosmology;
    ///
    /// let c = Cosmology::new(0.674, 0.315, 0.049, 0.0, 9.14e-5, 0.965, 0.811, -1.0, 0.0, 2.725);
    /// assert!(c.is_ok());
    /// ```
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        h: f64,
        omega_m: f64,
        omega_b: f64,
        omega_k: f64,
        omega_r: f64,
        n_s: f64,
        sigma_8: f64,
        w0: f64,
        wa: f64,
        t_cmb: f64,
    ) -> Result<Self, BrahmandaError> {
        if h <= 0.0 || !h.is_finite() {
            return Err(BrahmandaError::InvalidStructure(
                "Cosmology: h must be positive and finite".to_string(),
            ));
        }
        if omega_m < 0.0 || !omega_m.is_finite() {
            return Err(BrahmandaError::InvalidStructure(
                "Cosmology: omega_m must be non-negative and finite".to_string(),
            ));
        }
        if omega_b < 0.0 || omega_b > omega_m || !omega_b.is_finite() {
            return Err(BrahmandaError::InvalidStructure(
                "Cosmology: omega_b must be in [0, omega_m]".to_string(),
            ));
        }
        if sigma_8 <= 0.0 || !sigma_8.is_finite() {
            return Err(BrahmandaError::InvalidStructure(
                "Cosmology: sigma_8 must be positive and finite".to_string(),
            ));
        }
        if t_cmb <= 0.0 || !t_cmb.is_finite() {
            return Err(BrahmandaError::InvalidStructure(
                "Cosmology: t_cmb must be positive and finite".to_string(),
            ));
        }

        let h0_si = h * 100.0 / MPC_KM;
        let rho_crit = 3.0 * h0_si * h0_si / (8.0 * PI * G);
        let omega_lambda = 1.0 - omega_m - omega_k - omega_r;

        Ok(Self {
            h,
            omega_m,
            omega_b,
            omega_k,
            omega_r,
            n_s,
            sigma_8,
            w0,
            wa,
            t_cmb,
            h0_si,
            rho_crit,
            omega_lambda,
            a_s_cache: OnceLock::new(),
        })
    }

    /// Planck 2018 TT,TE,EE+lowE+lensing best-fit cosmology.
    ///
    /// ```
    /// use brahmanda::cosmology::Cosmology;
    ///
    /// let c = Cosmology::planck2018();
    /// assert!((c.sigma_8 - 0.811).abs() < 1e-10);
    /// ```
    #[must_use]
    pub fn planck2018() -> Self {
        Self::new(
            0.674, 0.315, 0.049, 0.0, 9.14e-5, 0.965, 0.811, -1.0, 0.0, 2.725,
        )
        .expect("Planck 2018 parameters are valid")
    }

    /// WMAP 9-year best-fit cosmology.
    ///
    /// ```
    /// use brahmanda::cosmology::Cosmology;
    ///
    /// let c = Cosmology::wmap9();
    /// assert!((c.h - 0.693).abs() < 1e-10);
    /// ```
    #[must_use]
    pub fn wmap9() -> Self {
        Self::new(
            0.693, 0.286, 0.047, 0.0, 8.5e-5, 0.965, 0.82, -1.0, 0.0, 2.725,
        )
        .expect("WMAP9 parameters are valid")
    }

    /// Power spectrum amplitude A_s, normalized so that σ(8 Mpc/h) = σ_8.
    ///
    /// Lazily computed and cached. The first call performs a numerical
    /// integration; subsequent calls return the cached value.
    ///
    /// ```
    /// use brahmanda::cosmology::Cosmology;
    ///
    /// let c = Cosmology::planck2018();
    /// let a_s = c.power_spectrum_amplitude();
    /// assert!(a_s > 0.0);
    /// ```
    pub fn power_spectrum_amplitude(&self) -> f64 {
        *self.a_s_cache.get_or_init(|| {
            // Compute σ8_raw² = ∫ k² k^n_s T²(k) W²(8k) dk / (2π²) in log-k space
            // u = ln(k), du = dk/k, so dk = k du
            // integrand = k³ k^n_s T²(k) W²(8k) / (2π²)
            let ln_k_min = -4.0 * std::f64::consts::LN_10; // ln(1e-4)
            let ln_k_max = 100.0_f64.ln();

            let integrand = |u: f64| -> f64 {
                let k = u.exp();
                let t = crate::power_spectrum::transfer_function_eh(
                    k,
                    self.omega_m,
                    self.omega_b,
                    self.h,
                )
                .unwrap_or(0.0);
                let w = FilterFunction::TopHat.window(k * 8.0);
                k.powi(3) * k.powf(self.n_s) * t * t * w * w / (2.0 * PI * PI)
            };

            let sigma8_raw_sq =
                hisab::calc::integral_adaptive_simpson(integrand, ln_k_min, ln_k_max, 1e-8)
                    .expect("sigma8 normalization integral must converge");

            if sigma8_raw_sq > 0.0 {
                self.sigma_8 * self.sigma_8 / sigma8_raw_sq
            } else {
                1.0
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_planck2018_preset() {
        let c = Cosmology::planck2018();
        assert!((c.h - 0.674).abs() < 1e-10);
        assert!((c.omega_m - 0.315).abs() < 1e-10);
        assert!((c.omega_b - 0.049).abs() < 1e-10);
        assert!((c.omega_lambda - 0.6849086).abs() < 1e-4);
        assert!(c.h0_si > 0.0);
        assert!(c.rho_crit > 0.0);
    }

    #[test]
    fn test_wmap9_preset() {
        let c = Cosmology::wmap9();
        assert!((c.h - 0.693).abs() < 1e-10);
        assert!((c.omega_m - 0.286).abs() < 1e-10);
    }

    #[test]
    fn test_invalid_h() {
        assert!(Cosmology::new(0.0, 0.3, 0.04, 0.0, 0.0, 0.96, 0.8, -1.0, 0.0, 2.7).is_err());
        assert!(Cosmology::new(-1.0, 0.3, 0.04, 0.0, 0.0, 0.96, 0.8, -1.0, 0.0, 2.7).is_err());
        assert!(Cosmology::new(f64::NAN, 0.3, 0.04, 0.0, 0.0, 0.96, 0.8, -1.0, 0.0, 2.7).is_err());
    }

    #[test]
    fn test_invalid_omega_b() {
        // omega_b > omega_m
        assert!(Cosmology::new(0.7, 0.3, 0.5, 0.0, 0.0, 0.96, 0.8, -1.0, 0.0, 2.7).is_err());
    }

    #[test]
    fn test_invalid_sigma_8() {
        assert!(Cosmology::new(0.7, 0.3, 0.04, 0.0, 0.0, 0.96, 0.0, -1.0, 0.0, 2.7).is_err());
        assert!(Cosmology::new(0.7, 0.3, 0.04, 0.0, 0.0, 0.96, -1.0, -1.0, 0.0, 2.7).is_err());
    }

    #[test]
    fn test_power_spectrum_amplitude_positive() {
        let c = Cosmology::planck2018();
        let a_s = c.power_spectrum_amplitude();
        assert!(a_s > 0.0, "A_s = {a_s}");
    }

    #[test]
    fn test_filter_tophat_at_zero() {
        assert!((FilterFunction::TopHat.window(0.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_filter_gaussian_at_zero() {
        assert!((FilterFunction::Gaussian.window(0.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_filter_sharpk() {
        assert!((FilterFunction::SharpK.window(0.5) - 1.0).abs() < 1e-10);
        assert!((FilterFunction::SharpK.window(1.5)).abs() < 1e-10);
    }

    #[test]
    fn test_clone_fresh_cache() {
        let c = Cosmology::planck2018();
        let _ = c.power_spectrum_amplitude(); // populate cache
        let c2 = c.clone();
        // The clone's cache is fresh (empty), but calling it should give the same value
        let a1 = c.power_spectrum_amplitude();
        let a2 = c2.power_spectrum_amplitude();
        assert!((a1 - a2).abs() / a1 < 1e-6);
    }

    #[test]
    fn test_serde_roundtrip() {
        let c = Cosmology::planck2018();
        let json = serde_json::to_string(&c).unwrap();
        let back: Cosmology = serde_json::from_str(&json).unwrap();
        assert_eq!(c, back);
    }
}
