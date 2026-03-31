//! Dark matter halos — density profiles, mass functions, concentration.
//!
//! NFW profile, halo mass function (Press-Schechter), virial relations.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::constants::{G, M_SUN, RHO_CRIT};
use crate::error::{ensure_finite, require_finite, BrahmandaError};

/// NFW (Navarro-Frenk-White) dark matter halo density profile.
///
/// ρ(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]
///
/// # Arguments
/// * `r` — Radius from halo center (kpc).
/// * `rho_s` — Characteristic density (M_sun/kpc³).
/// * `r_s` — Scale radius (kpc).
#[inline]
pub fn nfw_density(r: f64, rho_s: f64, r_s: f64) -> Result<f64, BrahmandaError> {
    require_finite(r, "nfw_density")?;
    require_finite(rho_s, "nfw_density")?;
    require_finite(r_s, "nfw_density")?;
    if r <= 0.0 || rho_s <= 0.0 || r_s <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "nfw_density: r, rho_s, r_s must be positive".to_string(),
        ));
    }
    let x = r / r_s;
    ensure_finite(rho_s / (x * (1.0 + x) * (1.0 + x)), "nfw_density")
}

/// NFW enclosed mass within radius r.
///
/// M(r) = 4π ρ_s r_s³ [ln(1 + r/r_s) - (r/r_s)/(1 + r/r_s)]
#[inline]
pub fn nfw_enclosed_mass(r: f64, rho_s: f64, r_s: f64) -> Result<f64, BrahmandaError> {
    require_finite(r, "nfw_enclosed_mass")?;
    require_finite(rho_s, "nfw_enclosed_mass")?;
    require_finite(r_s, "nfw_enclosed_mass")?;
    if r <= 0.0 || rho_s <= 0.0 || r_s <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "nfw_enclosed_mass: r, rho_s, r_s must be positive".to_string(),
        ));
    }
    let x = r / r_s;
    let factor = (1.0 + x).ln() - x / (1.0 + x);
    ensure_finite(
        4.0 * PI * rho_s * r_s * r_s * r_s * factor,
        "nfw_enclosed_mass",
    )
}

/// NFW circular velocity at radius r (km/s).
///
/// v_c(r) = √(GM(r)/r)
///
/// # Arguments
/// * `r_kpc` — Radius in kpc.
/// * `rho_s_msun_kpc3` — Characteristic density in M_sun/kpc³.
/// * `r_s_kpc` — Scale radius in kpc.
pub fn nfw_circular_velocity(
    r_kpc: f64,
    rho_s_msun_kpc3: f64,
    r_s_kpc: f64,
) -> Result<f64, BrahmandaError> {
    let m = nfw_enclosed_mass(r_kpc, rho_s_msun_kpc3, r_s_kpc)?;
    // Convert: M in M_sun, r in kpc → v in km/s
    let m_kg = m * M_SUN;
    let r_m = r_kpc * 3.085_677_581e19; // kpc → m
    let v_ms = (G * m_kg / r_m).sqrt();
    ensure_finite(v_ms / 1e3, "nfw_circular_velocity") // m/s → km/s
}

/// Virial radius for a halo of mass M_vir (kpc).
///
/// R_vir = [3 M_vir / (4π × Δ_vir × ρ_crit)]^(1/3)
///
/// Uses Δ_vir = 200 (standard definition).
pub fn virial_radius(m_vir_msun: f64) -> Result<f64, BrahmandaError> {
    require_finite(m_vir_msun, "virial_radius")?;
    if m_vir_msun <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "virial_radius: mass must be positive".to_string(),
        ));
    }
    let delta = 200.0;
    let m_kg = m_vir_msun * M_SUN;
    let r_m = (3.0 * m_kg / (4.0 * PI * delta * RHO_CRIT)).cbrt();
    let r_kpc = r_m / 3.085_677_581e19;
    ensure_finite(r_kpc, "virial_radius")
}

/// Concentration parameter: c = R_vir / r_s.
///
/// Typical: c ≈ 5–15 for halos of 10¹¹–10¹⁵ M_sun.
/// Uses Dutton & Maccio (2014) fitting formula:
/// log10(c) = 0.905 - 0.101 × log10(M_vir / (10¹² h⁻¹ M_sun))
pub fn concentration_dutton_maccio(m_vir_msun: f64) -> Result<f64, BrahmandaError> {
    require_finite(m_vir_msun, "concentration_dutton_maccio")?;
    if m_vir_msun <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "concentration_dutton_maccio: mass must be positive".to_string(),
        ));
    }
    let h = 0.674; // Planck 2018
    let log_m = (m_vir_msun / (1e12 / h)).log10();
    let log_c = 0.905 - 0.101 * log_m;
    ensure_finite(10.0_f64.powf(log_c), "concentration_dutton_maccio")
}

/// Halo properties bundle.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct HaloProperties {
    /// Virial mass (solar masses).
    pub m_vir_msun: f64,
    /// Virial radius (kpc).
    pub r_vir_kpc: f64,
    /// Scale radius (kpc).
    pub r_s_kpc: f64,
    /// Concentration c = R_vir / r_s.
    pub concentration: f64,
}

impl HaloProperties {
    /// Compute halo properties from virial mass.
    pub fn from_mass(m_vir_msun: f64) -> Result<Self, BrahmandaError> {
        require_finite(m_vir_msun, "HaloProperties::from_mass")?;
        let r_vir = virial_radius(m_vir_msun)?;
        let c = concentration_dutton_maccio(m_vir_msun)?;
        let r_s = r_vir / c;
        Ok(Self {
            m_vir_msun,
            r_vir_kpc: r_vir,
            r_s_kpc: r_s,
            concentration: c,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nfw_density_decreases() {
        let d1 = nfw_density(10.0, 1e7, 20.0).unwrap();
        let d2 = nfw_density(100.0, 1e7, 20.0).unwrap();
        assert!(d1 > d2, "NFW density should decrease outward");
    }

    #[test]
    fn test_nfw_mass_increases() {
        let m1 = nfw_enclosed_mass(10.0, 1e7, 20.0).unwrap();
        let m2 = nfw_enclosed_mass(100.0, 1e7, 20.0).unwrap();
        assert!(m2 > m1, "enclosed mass should increase with radius");
    }

    #[test]
    fn test_virial_radius_milky_way() {
        // MW: ~10¹² M_sun → R_vir ~ 200-300 kpc
        let r = virial_radius(1e12).unwrap();
        assert!(r > 150.0 && r < 350.0, "MW virial radius: {r} kpc");
    }

    #[test]
    fn test_concentration_milky_way() {
        let c = concentration_dutton_maccio(1e12).unwrap();
        assert!(c > 5.0 && c < 15.0, "MW concentration: {c}");
    }

    #[test]
    fn test_halo_properties() {
        let halo = HaloProperties::from_mass(1e12).unwrap();
        assert!(halo.r_vir_kpc > 0.0);
        assert!(halo.r_s_kpc > 0.0);
        assert!(halo.concentration > 0.0);
        assert!((halo.r_vir_kpc / halo.r_s_kpc - halo.concentration).abs() < 1e-6);
    }

    #[test]
    fn test_invalid_rejected() {
        assert!(nfw_density(-1.0, 1e7, 20.0).is_err());
        assert!(virial_radius(0.0).is_err());
        assert!(concentration_dutton_maccio(f64::NAN).is_err());
    }
}
