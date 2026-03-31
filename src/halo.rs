//! Dark matter halos — density profiles, mass functions, concentration.
//!
//! NFW profile, halo mass function (Press-Schechter), virial relations.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::constants::{G, KPC_M, M_SUN, OMEGA_M, RHO_CRIT};
use crate::error::{ensure_finite, require_finite, BrahmandaError};
use crate::power_spectrum;

/// NFW (Navarro-Frenk-White) dark matter halo density profile.
///
/// ρ(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]
///
/// # Arguments
/// * `r` — Radius from halo center (kpc).
/// * `rho_s` — Characteristic density (M_sun/kpc³).
/// * `r_s` — Scale radius (kpc).
///
/// ```
/// use brahmanda::halo::nfw_density;
///
/// let d_inner = nfw_density(10.0, 1e7, 20.0).unwrap();
/// let d_outer = nfw_density(100.0, 1e7, 20.0).unwrap();
/// assert!(d_inner > d_outer); // density decreases outward
/// ```
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
///
/// ```
/// use brahmanda::halo::nfw_enclosed_mass;
///
/// let m1 = nfw_enclosed_mass(10.0, 1e7, 20.0).unwrap();
/// let m2 = nfw_enclosed_mass(100.0, 1e7, 20.0).unwrap();
/// assert!(m2 > m1); // enclosed mass increases with radius
/// ```
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
///
/// ```
/// use brahmanda::halo::nfw_circular_velocity;
///
/// let v = nfw_circular_velocity(8.0, 1e7, 30.0).unwrap();
/// assert!(v > 0.0); // positive velocity
/// ```
pub fn nfw_circular_velocity(
    r_kpc: f64,
    rho_s_msun_kpc3: f64,
    r_s_kpc: f64,
) -> Result<f64, BrahmandaError> {
    let m = nfw_enclosed_mass(r_kpc, rho_s_msun_kpc3, r_s_kpc)?;
    // Convert: M in M_sun, r in kpc → v in km/s
    let m_kg = m * M_SUN;
    let r_m = r_kpc * KPC_M;
    let v_ms = (G * m_kg / r_m).sqrt();
    ensure_finite(v_ms / 1e3, "nfw_circular_velocity") // m/s → km/s
}

/// Virial radius for a halo of mass M_vir (kpc).
///
/// R_vir = [3 M_vir / (4π × Δ_vir × ρ_crit)]^(1/3)
///
/// Uses Δ_vir = 200 (standard definition).
///
/// ```
/// use brahmanda::halo::virial_radius;
///
/// // Milky Way-mass halo: ~200-300 kpc
/// let r = virial_radius(1e12).unwrap();
/// assert!(r > 150.0 && r < 350.0);
/// ```
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
    let r_kpc = r_m / KPC_M;
    ensure_finite(r_kpc, "virial_radius")
}

/// Concentration parameter: c = R_vir / r_s.
///
/// Typical: c ≈ 5–15 for halos of 10¹¹–10¹⁵ M_sun.
/// Uses Dutton & Maccio (2014) fitting formula:
/// log10(c) = 0.905 - 0.101 × log10(M_vir / (10¹² h⁻¹ M_sun))
///
/// ```
/// use brahmanda::halo::concentration_dutton_maccio;
///
/// let c = concentration_dutton_maccio(1e12).unwrap();
/// assert!(c > 5.0 && c < 15.0);
/// ```
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
#[must_use]
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
    ///
    /// ```
    /// use brahmanda::halo::HaloProperties;
    ///
    /// let halo = HaloProperties::from_mass(1e12).unwrap();
    /// assert!(halo.r_vir_kpc > 150.0 && halo.r_vir_kpc < 350.0);
    /// assert!(halo.concentration > 5.0 && halo.concentration < 15.0);
    /// ```
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

/// Critical overdensity for spherical collapse (EdS).
pub const DELTA_C: f64 = 1.686;

/// Lagrangian radius for a halo of mass M (Mpc/h).
///
/// R = (3M / 4π ρ̄_m)^(1/3), where ρ̄_m = Ω_m × ρ_crit.
///
/// ```
/// use brahmanda::halo::lagrangian_radius;
///
/// let r = lagrangian_radius(1e12).unwrap();
/// assert!(r > 0.0);
/// ```
pub fn lagrangian_radius(m_msun: f64) -> Result<f64, BrahmandaError> {
    require_finite(m_msun, "lagrangian_radius")?;
    if m_msun <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "lagrangian_radius: mass must be positive".to_string(),
        ));
    }
    let m_kg = m_msun * M_SUN;
    let rho_m = OMEGA_M * RHO_CRIT; // mean matter density (kg/m³)
    let r_m = (3.0 * m_kg / (4.0 * PI * rho_m)).cbrt();
    // Convert m → Mpc/h (h = 0.674)
    let h = 0.674;
    let mpc_m = 3.085_677_581e22;
    ensure_finite(r_m / mpc_m * h, "lagrangian_radius")
}

/// Peak height ν = δ_c / σ(M, z).
///
/// The dimensionless threshold that determines how rare a halo of mass M
/// is at redshift z. Higher ν means rarer (more massive) halos.
///
/// ```
/// use brahmanda::halo::peak_height;
///
/// let nu = peak_height(1e12, 0.0).unwrap();
/// assert!(nu > 0.0);
/// ```
pub fn peak_height(m_msun: f64, z: f64) -> Result<f64, BrahmandaError> {
    let r = lagrangian_radius(m_msun)?;
    let sigma = power_spectrum::sigma_r(r, z)?;
    ensure_finite(DELTA_C / sigma, "peak_height")
}

/// Press-Schechter halo mass function: dn/dlnM.
///
/// Returns the comoving number density of halos per unit ln(M)
/// in units of h³/Mpc³.
///
/// dn/dlnM = (ρ̄_m / M) × |dlnσ/dlnM| × f(ν)
///
/// where f(ν) = √(2/π) × ν × exp(-ν²/2) is the PS multiplicity function.
///
/// # Arguments
/// * `m_msun` — Halo mass in solar masses.
/// * `z` — Redshift.
///
/// ```
/// use brahmanda::halo::press_schechter_dndlnm;
///
/// // More low-mass halos than high-mass
/// let n_low = press_schechter_dndlnm(1e10, 0.0).unwrap();
/// let n_high = press_schechter_dndlnm(1e14, 0.0).unwrap();
/// assert!(n_low > n_high);
/// ```
pub fn press_schechter_dndlnm(m_msun: f64, z: f64) -> Result<f64, BrahmandaError> {
    require_finite(m_msun, "press_schechter_dndlnm")?;
    require_finite(z, "press_schechter_dndlnm")?;
    if m_msun <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "press_schechter_dndlnm: mass must be positive".to_string(),
        ));
    }

    let nu = peak_height(m_msun, z)?;

    // PS multiplicity: f(ν) = √(2/π) × ν × exp(-ν²/2)
    let f_nu = (2.0 / PI).sqrt() * nu * (-nu * nu / 2.0).exp();

    // |dlnσ/dlnM| ≈ (n_s + 3) / 6 for our power-law σ(R) approximation
    // Since σ ∝ R^(-(n_s+3)/6) and M ∝ R³, dlnσ/dlnM = -(n_s+3)/18
    let n_s = crate::constants::N_S;
    let dlnsigma_dlnm = (n_s + 3.0) / 18.0; // absolute value

    // Mean matter density in h²M_sun/Mpc³
    let h = 0.674;
    let rho_m_kg_m3 = OMEGA_M * RHO_CRIT;
    let mpc_m = 3.085_677_581e22;
    let rho_m_msun_mpc3 = rho_m_kg_m3 * (mpc_m * mpc_m * mpc_m) / M_SUN;
    let rho_m_h2 = rho_m_msun_mpc3 / (h * h * h); // h³ Mpc⁻³ M_sun

    let dndlnm = (rho_m_h2 / m_msun) * dlnsigma_dlnm * f_nu;
    ensure_finite(dndlnm, "press_schechter_dndlnm")
}

/// Linear halo bias — Mo & White (1996).
///
/// b(M, z) = 1 + (ν² - 1) / δ_c
///
/// where ν = δ_c / σ(M, z) is the peak height.
///
/// Massive halos (ν > 1) are positively biased (b > 1),
/// while low-mass halos can be anti-biased (b < 1).
///
/// ```
/// use brahmanda::halo::bias_mo_white;
///
/// // Cluster-mass halos are strongly biased
/// let b_cluster = bias_mo_white(1e15, 0.0).unwrap();
/// assert!(b_cluster > 2.0);
///
/// // MW-mass halos are mildly biased
/// let b_mw = bias_mo_white(1e12, 0.0).unwrap();
/// assert!(b_mw > 0.5 && b_mw < 3.0);
/// ```
pub fn bias_mo_white(m_msun: f64, z: f64) -> Result<f64, BrahmandaError> {
    let nu = peak_height(m_msun, z)?;
    ensure_finite(1.0 + (nu * nu - 1.0) / DELTA_C, "bias_mo_white")
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

    #[test]
    fn test_lagrangian_radius_increases_with_mass() {
        let r1 = lagrangian_radius(1e10).unwrap();
        let r2 = lagrangian_radius(1e12).unwrap();
        let r3 = lagrangian_radius(1e14).unwrap();
        assert!(r3 > r2 && r2 > r1);
    }

    #[test]
    fn test_peak_height_increases_with_mass() {
        // More massive halos → higher ν (rarer peaks)
        let nu1 = peak_height(1e10, 0.0).unwrap();
        let nu2 = peak_height(1e14, 0.0).unwrap();
        assert!(nu2 > nu1, "ν should increase with mass");
    }

    #[test]
    fn test_press_schechter_decreases_with_mass() {
        let n1 = press_schechter_dndlnm(1e10, 0.0).unwrap();
        let n2 = press_schechter_dndlnm(1e12, 0.0).unwrap();
        let n3 = press_schechter_dndlnm(1e14, 0.0).unwrap();
        assert!(n1 > n2, "more low-mass halos");
        assert!(n2 > n3, "fewer high-mass halos");
    }

    #[test]
    fn test_press_schechter_positive() {
        for m in [1e8, 1e10, 1e12, 1e14] {
            let n = press_schechter_dndlnm(m, 0.0).unwrap();
            assert!(n > 0.0, "dn/dlnM should be positive for M={m}");
        }
    }

    #[test]
    fn test_bias_cluster_strongly_biased() {
        let b = bias_mo_white(1e15, 0.0).unwrap();
        assert!(b > 2.0, "cluster bias: {b}");
    }

    #[test]
    fn test_bias_increases_with_mass() {
        let b1 = bias_mo_white(1e11, 0.0).unwrap();
        let b2 = bias_mo_white(1e13, 0.0).unwrap();
        let b3 = bias_mo_white(1e15, 0.0).unwrap();
        assert!(b3 > b2 && b2 > b1, "bias should increase with mass");
    }
}
