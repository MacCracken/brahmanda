//! Dark matter halos — density profiles, mass functions, concentration, accretion.
//!
//! NFW profile, halo mass functions (Press-Schechter, Sheth-Tormen, Tinker08),
//! virial relations, mass accretion history, subhalo abundance matching.
//!
//! All functions that depend on cosmological parameters take a
//! [`Cosmology`] reference.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::constants::{G, KPC_M, M_SUN, MPC_M};
use crate::cosmology::Cosmology;
use crate::error::{BrahmandaError, ensure_finite, require_finite};
use crate::power_spectrum;

/// Mean matter density in M_sun / (Mpc/h)³.
fn mean_matter_density_h2(cosmo: &Cosmology) -> f64 {
    let rho_m_kg_m3 = cosmo.omega_m * cosmo.rho_crit;
    let rho_m_msun_mpc3 = rho_m_kg_m3 * (MPC_M * MPC_M * MPC_M) / M_SUN;
    rho_m_msun_mpc3 / (cosmo.h * cosmo.h * cosmo.h)
}

/// NFW (Navarro-Frenk-White) dark matter halo density profile.
///
/// ρ(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]
///
/// ```
/// use brahmanda::halo::nfw_density;
///
/// let d_inner = nfw_density(10.0, 1e7, 20.0).unwrap();
/// let d_outer = nfw_density(100.0, 1e7, 20.0).unwrap();
/// assert!(d_inner > d_outer);
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
/// assert!(m2 > m1);
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
/// ```
/// use brahmanda::halo::nfw_circular_velocity;
///
/// let v = nfw_circular_velocity(8.0, 1e7, 30.0).unwrap();
/// assert!(v > 0.0);
/// ```
pub fn nfw_circular_velocity(
    r_kpc: f64,
    rho_s_msun_kpc3: f64,
    r_s_kpc: f64,
) -> Result<f64, BrahmandaError> {
    let m = nfw_enclosed_mass(r_kpc, rho_s_msun_kpc3, r_s_kpc)?;
    let m_kg = m * M_SUN;
    let r_m = r_kpc * KPC_M;
    let v_ms = (G * m_kg / r_m).sqrt();
    ensure_finite(v_ms / 1e3, "nfw_circular_velocity")
}

/// Virial radius for a halo of mass M_vir (kpc).
///
/// R_vir = [3 M_vir / (4π × Δ_vir × ρ_crit)]^(1/3)
///
/// Uses Δ_vir = 200 (standard definition).
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::halo::virial_radius;
///
/// let cosmo = Cosmology::planck2018();
/// let r = virial_radius(1e12, &cosmo).unwrap();
/// assert!(r > 150.0 && r < 350.0);
/// ```
pub fn virial_radius(m_vir_msun: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(m_vir_msun, "virial_radius")?;
    if m_vir_msun <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "virial_radius: mass must be positive".to_string(),
        ));
    }
    let delta = 200.0;
    let m_kg = m_vir_msun * M_SUN;
    let r_m = (3.0 * m_kg / (4.0 * PI * delta * cosmo.rho_crit)).cbrt();
    let r_kpc = r_m / KPC_M;
    ensure_finite(r_kpc, "virial_radius")
}

/// Concentration parameter: c = R_vir / r_s.
///
/// Uses Dutton & Maccio (2014) fitting formula.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::halo::concentration_dutton_maccio;
///
/// let cosmo = Cosmology::planck2018();
/// let c = concentration_dutton_maccio(1e12, &cosmo).unwrap();
/// assert!(c > 5.0 && c < 15.0);
/// ```
pub fn concentration_dutton_maccio(
    m_vir_msun: f64,
    cosmo: &Cosmology,
) -> Result<f64, BrahmandaError> {
    require_finite(m_vir_msun, "concentration_dutton_maccio")?;
    if m_vir_msun <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "concentration_dutton_maccio: mass must be positive".to_string(),
        ));
    }
    let log_m = (m_vir_msun / (1e12 / cosmo.h)).log10();
    let log_c = 0.905 - 0.101 * log_m;
    ensure_finite(10.0_f64.powf(log_c), "concentration_dutton_maccio")
}

/// Halo properties bundle.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
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
    /// use brahmanda::Cosmology;
    /// use brahmanda::halo::HaloProperties;
    ///
    /// let cosmo = Cosmology::planck2018();
    /// let halo = HaloProperties::from_mass(1e12, &cosmo).unwrap();
    /// assert!(halo.r_vir_kpc > 150.0 && halo.r_vir_kpc < 350.0);
    /// assert!(halo.concentration > 5.0 && halo.concentration < 15.0);
    /// ```
    pub fn from_mass(m_vir_msun: f64, cosmo: &Cosmology) -> Result<Self, BrahmandaError> {
        require_finite(m_vir_msun, "HaloProperties::from_mass")?;
        let r_vir = virial_radius(m_vir_msun, cosmo)?;
        let c = concentration_dutton_maccio(m_vir_msun, cosmo)?;
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
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::halo::lagrangian_radius;
///
/// let cosmo = Cosmology::planck2018();
/// let r = lagrangian_radius(1e12, &cosmo).unwrap();
/// assert!(r > 0.0);
/// ```
pub fn lagrangian_radius(m_msun: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(m_msun, "lagrangian_radius")?;
    if m_msun <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "lagrangian_radius: mass must be positive".to_string(),
        ));
    }
    let m_kg = m_msun * M_SUN;
    let rho_m = cosmo.omega_m * cosmo.rho_crit;
    let r_m = (3.0 * m_kg / (4.0 * PI * rho_m)).cbrt();
    ensure_finite(r_m / MPC_M * cosmo.h, "lagrangian_radius")
}

/// Peak height ν = δ_c / σ(M, z).
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::halo::peak_height;
///
/// let cosmo = Cosmology::planck2018();
/// let nu = peak_height(1e12, 0.0, &cosmo).unwrap();
/// assert!(nu > 0.0);
/// ```
pub fn peak_height(m_msun: f64, z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    let r = lagrangian_radius(m_msun, cosmo)?;
    let sigma = power_spectrum::sigma_r(r, z, cosmo)?;
    ensure_finite(DELTA_C / sigma, "peak_height")
}

/// Convert a multiplicity function f(ν) value to dn/dlnM (h³/Mpc³).
///
/// dn/dlnM = (ρ̄_m / M) × |dlnσ/dlnM| × f(ν)
///
/// Computes exact |dlnσ/dlnM| via finite-difference of σ(R).
fn mass_function_from_fnu(
    f_nu: f64,
    m_msun: f64,
    z: f64,
    cosmo: &Cosmology,
) -> Result<f64, BrahmandaError> {
    let r = lagrangian_radius(m_msun, cosmo)?;
    let dr = r * 0.01;
    let s_lo = power_spectrum::sigma_r(r - dr, z, cosmo)?;
    let s_hi = power_spectrum::sigma_r(r + dr, z, cosmo)?;
    let dlns_dlnr = (s_hi.ln() - s_lo.ln()) / ((r + dr).ln() - (r - dr).ln());
    // dlnσ/dlnM = dlnσ/dlnR × dlnR/dlnM = dlnσ/dlnR / 3
    let dlns_dlnm = dlns_dlnr / 3.0;

    let rho_m = mean_matter_density_h2(cosmo);
    ensure_finite(
        (rho_m / m_msun) * dlns_dlnm.abs() * f_nu,
        "mass_function_from_fnu",
    )
}

/// Press-Schechter halo mass function: dn/dlnM.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::halo::press_schechter_dndlnm;
///
/// let cosmo = Cosmology::planck2018();
/// let n_low = press_schechter_dndlnm(1e10, 0.0, &cosmo).unwrap();
/// let n_high = press_schechter_dndlnm(1e14, 0.0, &cosmo).unwrap();
/// assert!(n_low > n_high);
/// ```
pub fn press_schechter_dndlnm(
    m_msun: f64,
    z: f64,
    cosmo: &Cosmology,
) -> Result<f64, BrahmandaError> {
    require_finite(m_msun, "press_schechter_dndlnm")?;
    require_finite(z, "press_schechter_dndlnm")?;
    if m_msun <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "press_schechter_dndlnm: mass must be positive".to_string(),
        ));
    }

    let nu = peak_height(m_msun, z, cosmo)?;
    let f_nu = (2.0 / PI).sqrt() * nu * (-nu * nu / 2.0).exp();
    mass_function_from_fnu(f_nu, m_msun, z, cosmo)
}

/// Linear halo bias — Mo & White (1996).
///
/// b(M, z) = 1 + (ν² - 1) / δ_c
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::halo::bias_mo_white;
///
/// let cosmo = Cosmology::planck2018();
/// let b = bias_mo_white(1e15, 0.0, &cosmo).unwrap();
/// assert!(b > 2.0);
/// ```
pub fn bias_mo_white(m_msun: f64, z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    let nu = peak_height(m_msun, z, cosmo)?;
    ensure_finite(1.0 + (nu * nu - 1.0) / DELTA_C, "bias_mo_white")
}

// ============================================================
// Sheth-Tormen mass function
// ============================================================

/// Sheth-Tormen (1999) ellipsoidal collapse parameters.
pub const ST_A: f64 = 0.707;
/// Sheth-Tormen p parameter.
pub const ST_P: f64 = 0.3;
/// Sheth-Tormen normalization.
pub const ST_BIG_A: f64 = 0.3222;

/// Sheth-Tormen halo mass function: dn/dlnM.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::halo::{sheth_tormen_dndlnm, press_schechter_dndlnm};
///
/// let cosmo = Cosmology::planck2018();
/// let n_st = sheth_tormen_dndlnm(1e15, 0.0, &cosmo).unwrap();
/// let n_ps = press_schechter_dndlnm(1e15, 0.0, &cosmo).unwrap();
/// assert!(n_st > n_ps);
/// ```
pub fn sheth_tormen_dndlnm(m_msun: f64, z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(m_msun, "sheth_tormen_dndlnm")?;
    require_finite(z, "sheth_tormen_dndlnm")?;
    if m_msun <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "sheth_tormen_dndlnm: mass must be positive".to_string(),
        ));
    }

    let nu = peak_height(m_msun, z, cosmo)?;
    let a_nu2 = ST_A * nu * nu;
    let f_nu =
        ST_BIG_A * (2.0 * ST_A / PI).sqrt() * nu * (1.0 + a_nu2.powf(-ST_P)) * (-a_nu2 / 2.0).exp();

    mass_function_from_fnu(f_nu, m_msun, z, cosmo)
}

/// Sheth-Tormen linear halo bias.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::halo::bias_sheth_tormen;
///
/// let cosmo = Cosmology::planck2018();
/// let b = bias_sheth_tormen(1e14, 0.0, &cosmo).unwrap();
/// assert!(b > 1.0);
/// ```
pub fn bias_sheth_tormen(m_msun: f64, z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    let nu = peak_height(m_msun, z, cosmo)?;
    let a_nu2 = ST_A * nu * nu;
    let b = 1.0 + (a_nu2 - 1.0) / DELTA_C + 2.0 * ST_P / (DELTA_C * (1.0 + a_nu2.powf(ST_P)));
    ensure_finite(b, "bias_sheth_tormen")
}

// ============================================================
// Tinker (2008) mass function & (2010) bias
// ============================================================

/// Tinker et al. (2008) halo mass function: dn/dlnM.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::halo::tinker08_dndlnm;
///
/// let cosmo = Cosmology::planck2018();
/// let n_low = tinker08_dndlnm(1e10, 0.0, &cosmo).unwrap();
/// let n_high = tinker08_dndlnm(1e14, 0.0, &cosmo).unwrap();
/// assert!(n_low > n_high);
/// ```
pub fn tinker08_dndlnm(m_msun: f64, z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(m_msun, "tinker08_dndlnm")?;
    require_finite(z, "tinker08_dndlnm")?;
    if m_msun <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "tinker08_dndlnm: mass must be positive".to_string(),
        ));
    }

    let nu = peak_height(m_msun, z, cosmo)?;
    let sigma = DELTA_C / nu;

    let a0 = 0.186;
    let a_a = 1.47;
    let a_b = 2.57;
    let a_c = 1.19;

    let big_a = a0 * (1.0 + z).powf(-0.14);
    let a = a_a * (1.0 + z).powf(-0.06);
    let log_alpha = -(0.75 / ((z / 2.3).powi(2) + 1.0).log10()).max(-10.0);
    let b = a_b * (1.0 + z).powf(-10.0_f64.powf(log_alpha));
    let c = a_c;

    let f_sigma = big_a * ((sigma / b).powf(-a) + 1.0) * (-c / (sigma * sigma)).exp();

    mass_function_from_fnu(f_sigma, m_msun, z, cosmo)
}

/// Tinker et al. (2010) linear halo bias.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::halo::bias_tinker10;
///
/// let cosmo = Cosmology::planck2018();
/// let b = bias_tinker10(1e14, 0.0, &cosmo).unwrap();
/// assert!(b > 1.0);
/// ```
pub fn bias_tinker10(m_msun: f64, z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    let nu = peak_height(m_msun, z, cosmo)?;

    let big_a = 1.005;
    let a_a = 0.378;
    let big_b = 0.183;
    let a_b = 1.5;
    let big_c = 0.265;
    let a_c = 2.4;

    let nu_a = nu.powf(a_a);
    let b = 1.0 - big_a * nu_a / (nu_a + DELTA_C.powf(a_a))
        + big_b * nu.powf(a_b)
        + big_c * nu.powf(a_c);
    ensure_finite(b, "bias_tinker10")
}

// ============================================================
// Halo mass accretion history
// ============================================================

/// Halo mass accretion history — Wechsler et al. (2002).
///
/// M(z) = M₀ × exp(-α × z)
///
/// ```
/// use brahmanda::halo::mass_accretion_wechsler;
///
/// let m0 = 1e12;
/// let m = mass_accretion_wechsler(m0, 0.0, 1.0).unwrap();
/// assert!((m - m0).abs() < 1.0);
/// ```
pub fn mass_accretion_wechsler(m0_msun: f64, z: f64, alpha: f64) -> Result<f64, BrahmandaError> {
    require_finite(m0_msun, "mass_accretion_wechsler")?;
    require_finite(z, "mass_accretion_wechsler")?;
    require_finite(alpha, "mass_accretion_wechsler")?;
    if m0_msun <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "mass_accretion_wechsler: mass must be positive".to_string(),
        ));
    }
    if alpha <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "mass_accretion_wechsler: alpha must be positive".to_string(),
        ));
    }
    ensure_finite(m0_msun * (-alpha * z).exp(), "mass_accretion_wechsler")
}

/// Accretion rate parameter α from concentration.
///
/// ```
/// use brahmanda::halo::accretion_rate_from_concentration;
///
/// let alpha = accretion_rate_from_concentration(8.0).unwrap();
/// assert!(alpha > 1.0 && alpha < 3.0);
/// ```
#[inline]
pub fn accretion_rate_from_concentration(c: f64) -> Result<f64, BrahmandaError> {
    require_finite(c, "accretion_rate_from_concentration")?;
    if c <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "accretion_rate_from_concentration: c must be positive".to_string(),
        ));
    }
    ensure_finite(c / 4.1, "accretion_rate_from_concentration")
}

/// Halo mass accretion history — McBride et al. (2009).
///
/// M(z) = M₀ × (1+z)^β × exp(-γ z)
///
/// ```
/// use brahmanda::halo::mass_accretion_mcbride;
///
/// let m0 = 1e12;
/// let m = mass_accretion_mcbride(m0, 0.0, 0.1, 0.8).unwrap();
/// assert!((m - m0).abs() < 1.0);
/// ```
pub fn mass_accretion_mcbride(
    m0_msun: f64,
    z: f64,
    beta: f64,
    gamma: f64,
) -> Result<f64, BrahmandaError> {
    require_finite(m0_msun, "mass_accretion_mcbride")?;
    require_finite(z, "mass_accretion_mcbride")?;
    require_finite(beta, "mass_accretion_mcbride")?;
    require_finite(gamma, "mass_accretion_mcbride")?;
    if m0_msun <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "mass_accretion_mcbride: mass must be positive".to_string(),
        ));
    }
    if gamma <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "mass_accretion_mcbride: gamma must be positive".to_string(),
        ));
    }
    ensure_finite(
        m0_msun * (1.0 + z).powf(beta) * (-gamma * z).exp(),
        "mass_accretion_mcbride",
    )
}

// ============================================================
// Subhalo abundance matching (SHAM)
// ============================================================

/// Schechter stellar mass function: φ(M*) dn/dlog₁₀M*.
///
/// ```
/// use brahmanda::halo::schechter_smf;
///
/// let phi = schechter_smf(10.0, 10.65, 4.5e-3, -1.2).unwrap();
/// assert!(phi > 0.0);
/// ```
pub fn schechter_smf(
    log_m_star: f64,
    log_m_star_0: f64,
    phi_star: f64,
    alpha: f64,
) -> Result<f64, BrahmandaError> {
    require_finite(log_m_star, "schechter_smf")?;
    require_finite(log_m_star_0, "schechter_smf")?;
    require_finite(phi_star, "schechter_smf")?;
    require_finite(alpha, "schechter_smf")?;
    if phi_star <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "schechter_smf: phi_star must be positive".to_string(),
        ));
    }

    let x = 10.0_f64.powf(log_m_star - log_m_star_0);
    let phi = (10.0_f64.ln()) * phi_star * x.powf(alpha + 1.0) * (-x).exp();
    ensure_finite(phi, "schechter_smf")
}

/// Cumulative stellar mass function n(>M*).
///
/// ```
/// use brahmanda::halo::schechter_cumulative;
///
/// let n1 = schechter_cumulative(10.0, 10.65, 4.5e-3, -1.2).unwrap();
/// let n2 = schechter_cumulative(11.0, 10.65, 4.5e-3, -1.2).unwrap();
/// assert!(n1 > n2);
/// ```
pub fn schechter_cumulative(
    log_m_star_min: f64,
    log_m_star_0: f64,
    phi_star: f64,
    alpha: f64,
) -> Result<f64, BrahmandaError> {
    require_finite(log_m_star_min, "schechter_cumulative")?;
    let log_max = log_m_star_0 + 4.0;
    if log_m_star_min >= log_max {
        return Ok(0.0);
    }

    let n = 200_usize;
    let dlog = (log_max - log_m_star_min) / n as f64;
    let mut sum = 0.0;

    for i in 0..n {
        let l0 = log_m_star_min + i as f64 * dlog;
        let l1 = l0 + dlog / 2.0;
        let l2 = l0 + dlog;
        let f0 = schechter_smf(l0, log_m_star_0, phi_star, alpha)?;
        let f1 = schechter_smf(l1, log_m_star_0, phi_star, alpha)?;
        let f2 = schechter_smf(l2, log_m_star_0, phi_star, alpha)?;
        sum += (dlog / 6.0) * (f0 + 4.0 * f1 + f2);
    }

    ensure_finite(sum, "schechter_cumulative")
}

/// Subhalo abundance matching — stellar mass from halo mass.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::halo::sham_stellar_mass;
///
/// let cosmo = Cosmology::planck2018();
/// let log_ms = sham_stellar_mass(1e12, 0.0, 10.65, 4.5e-3, -1.2, &cosmo).unwrap();
/// assert!(log_ms > 9.0 && log_ms < 12.0);
/// ```
pub fn sham_stellar_mass(
    m_halo_msun: f64,
    z: f64,
    log_m_star_0: f64,
    phi_star: f64,
    alpha: f64,
    cosmo: &Cosmology,
) -> Result<f64, BrahmandaError> {
    require_finite(m_halo_msun, "sham_stellar_mass")?;
    require_finite(z, "sham_stellar_mass")?;
    if m_halo_msun <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "sham_stellar_mass: halo mass must be positive".to_string(),
        ));
    }

    let log_m_h = m_halo_msun.log10();
    let log_m_max = 16.0;
    let n_steps = 100_usize;
    let dlog = (log_m_max - log_m_h) / n_steps as f64;
    let mut n_halo = 0.0;

    for i in 0..n_steps {
        let l0 = log_m_h + i as f64 * dlog;
        let l1 = l0 + dlog / 2.0;
        let l2 = l0 + dlog;
        let f = |log_m: f64| -> Result<f64, BrahmandaError> {
            let m = 10.0_f64.powf(log_m);
            let dn = press_schechter_dndlnm(m, z, cosmo)?;
            Ok(dn * 10.0_f64.ln())
        };
        let f0 = f(l0)?;
        let f1 = f(l1)?;
        let f2 = f(l2)?;
        n_halo += (dlog / 6.0) * (f0 + 4.0 * f1 + f2);
    }

    if n_halo <= 0.0 {
        return Err(BrahmandaError::Computation(
            "sham_stellar_mass: halo cumulative density is zero".to_string(),
        ));
    }

    let mut lo = 6.0_f64;
    let mut hi = log_m_star_0 + 3.0;

    for _ in 0..80 {
        let mid = (lo + hi) / 2.0;
        let n_gal = schechter_cumulative(mid, log_m_star_0, phi_star, alpha)?;
        if n_gal > n_halo {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    ensure_finite((lo + hi) / 2.0, "sham_stellar_mass")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Cosmology;

    fn cosmo() -> Cosmology {
        Cosmology::planck2018()
    }

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
        let c = cosmo();
        let r = virial_radius(1e12, &c).unwrap();
        assert!(r > 150.0 && r < 350.0, "MW virial radius: {r} kpc");
    }

    #[test]
    fn test_concentration_milky_way() {
        let c = cosmo();
        let conc = concentration_dutton_maccio(1e12, &c).unwrap();
        assert!(conc > 5.0 && conc < 15.0, "MW concentration: {conc}");
    }

    #[test]
    fn test_halo_properties() {
        let c = cosmo();
        let halo = HaloProperties::from_mass(1e12, &c).unwrap();
        assert!(halo.r_vir_kpc > 0.0);
        assert!(halo.r_s_kpc > 0.0);
        assert!(halo.concentration > 0.0);
        assert!((halo.r_vir_kpc / halo.r_s_kpc - halo.concentration).abs() < 1e-6);
    }

    #[test]
    fn test_invalid_rejected() {
        let c = cosmo();
        assert!(nfw_density(-1.0, 1e7, 20.0).is_err());
        assert!(virial_radius(0.0, &c).is_err());
        assert!(concentration_dutton_maccio(f64::NAN, &c).is_err());
    }

    #[test]
    fn test_lagrangian_radius_increases_with_mass() {
        let c = cosmo();
        let r1 = lagrangian_radius(1e10, &c).unwrap();
        let r2 = lagrangian_radius(1e12, &c).unwrap();
        let r3 = lagrangian_radius(1e14, &c).unwrap();
        assert!(r3 > r2 && r2 > r1);
    }

    #[test]
    fn test_peak_height_increases_with_mass() {
        let c = cosmo();
        let nu1 = peak_height(1e10, 0.0, &c).unwrap();
        let nu2 = peak_height(1e14, 0.0, &c).unwrap();
        assert!(nu2 > nu1, "ν should increase with mass");
    }

    #[test]
    fn test_press_schechter_decreases_with_mass() {
        let c = cosmo();
        let n1 = press_schechter_dndlnm(1e10, 0.0, &c).unwrap();
        let n2 = press_schechter_dndlnm(1e12, 0.0, &c).unwrap();
        let n3 = press_schechter_dndlnm(1e14, 0.0, &c).unwrap();
        assert!(n1 > n2, "more low-mass halos");
        assert!(n2 > n3, "fewer high-mass halos");
    }

    #[test]
    fn test_press_schechter_positive() {
        let c = cosmo();
        for m in [1e8, 1e10, 1e12, 1e14] {
            let n = press_schechter_dndlnm(m, 0.0, &c).unwrap();
            assert!(n > 0.0, "dn/dlnM should be positive for M={m}");
        }
    }

    #[test]
    fn test_bias_cluster_strongly_biased() {
        let c = cosmo();
        let b = bias_mo_white(1e15, 0.0, &c).unwrap();
        assert!(b > 2.0, "cluster bias: {b}");
    }

    #[test]
    fn test_bias_increases_with_mass() {
        let c = cosmo();
        let b1 = bias_mo_white(1e11, 0.0, &c).unwrap();
        let b2 = bias_mo_white(1e13, 0.0, &c).unwrap();
        let b3 = bias_mo_white(1e15, 0.0, &c).unwrap();
        assert!(b3 > b2 && b2 > b1, "bias should increase with mass");
    }

    // -- Sheth-Tormen --

    #[test]
    fn test_sheth_tormen_positive() {
        let c = cosmo();
        for m in [1e8, 1e10, 1e12, 1e14] {
            let n = sheth_tormen_dndlnm(m, 0.0, &c).unwrap();
            assert!(n > 0.0, "ST dn/dlnM should be positive for M={m}");
        }
    }

    #[test]
    fn test_sheth_tormen_decreases_with_mass() {
        let c = cosmo();
        let n1 = sheth_tormen_dndlnm(1e10, 0.0, &c).unwrap();
        let n2 = sheth_tormen_dndlnm(1e12, 0.0, &c).unwrap();
        let n3 = sheth_tormen_dndlnm(1e14, 0.0, &c).unwrap();
        assert!(n1 > n2 && n2 > n3);
    }

    #[test]
    fn test_sheth_tormen_more_massive_than_ps() {
        let c = cosmo();
        let st = sheth_tormen_dndlnm(1e15, 0.0, &c).unwrap();
        let ps = press_schechter_dndlnm(1e15, 0.0, &c).unwrap();
        assert!(st > ps, "ST={st} should exceed PS={ps} at high mass");
    }

    #[test]
    fn test_bias_sheth_tormen_increases_with_mass() {
        let c = cosmo();
        let b1 = bias_sheth_tormen(1e11, 0.0, &c).unwrap();
        let b2 = bias_sheth_tormen(1e13, 0.0, &c).unwrap();
        let b3 = bias_sheth_tormen(1e15, 0.0, &c).unwrap();
        assert!(b3 > b2 && b2 > b1, "ST bias should increase with mass");
    }

    // -- Tinker 2008 mass function --

    #[test]
    fn test_tinker08_positive() {
        let c = cosmo();
        for m in [1e8, 1e10, 1e12, 1e14] {
            let n = tinker08_dndlnm(m, 0.0, &c).unwrap();
            assert!(n > 0.0, "Tinker08 dn/dlnM should be positive for M={m}");
        }
    }

    #[test]
    fn test_tinker08_decreases_with_mass() {
        let c = cosmo();
        let n1 = tinker08_dndlnm(1e10, 0.0, &c).unwrap();
        let n2 = tinker08_dndlnm(1e12, 0.0, &c).unwrap();
        let n3 = tinker08_dndlnm(1e14, 0.0, &c).unwrap();
        assert!(n1 > n2 && n2 > n3);
    }

    #[test]
    fn test_tinker08_invalid() {
        let c = cosmo();
        assert!(tinker08_dndlnm(0.0, 0.0, &c).is_err());
        assert!(tinker08_dndlnm(-1e12, 0.0, &c).is_err());
        assert!(tinker08_dndlnm(f64::NAN, 0.0, &c).is_err());
    }

    // -- Tinker 2010 bias --

    #[test]
    fn test_bias_tinker10_increases_with_mass() {
        let c = cosmo();
        let b1 = bias_tinker10(1e11, 0.0, &c).unwrap();
        let b2 = bias_tinker10(1e13, 0.0, &c).unwrap();
        let b3 = bias_tinker10(1e15, 0.0, &c).unwrap();
        assert!(
            b3 > b2 && b2 > b1,
            "Tinker10 bias should increase with mass"
        );
    }

    #[test]
    fn test_bias_tinker10_cluster_biased() {
        let c = cosmo();
        let b = bias_tinker10(1e14, 0.0, &c).unwrap();
        assert!(b > 1.0, "cluster bias = {b}");
    }

    // -- Mass accretion history --

    #[test]
    fn test_wechsler_z0_equals_m0() {
        let m = mass_accretion_wechsler(1e12, 0.0, 1.0).unwrap();
        assert!((m - 1e12).abs() < 1.0);
    }

    #[test]
    fn test_wechsler_decreases_with_z() {
        let m0 = mass_accretion_wechsler(1e12, 0.0, 1.0).unwrap();
        let m1 = mass_accretion_wechsler(1e12, 1.0, 1.0).unwrap();
        let m3 = mass_accretion_wechsler(1e12, 3.0, 1.0).unwrap();
        assert!(m0 > m1 && m1 > m3);
    }

    #[test]
    fn test_wechsler_higher_alpha_faster_growth() {
        let m_slow = mass_accretion_wechsler(1e12, 2.0, 0.5).unwrap();
        let m_fast = mass_accretion_wechsler(1e12, 2.0, 1.5).unwrap();
        assert!(m_slow > m_fast, "higher α → lower M(z>0)");
    }

    #[test]
    fn test_accretion_rate_from_concentration() {
        let alpha = accretion_rate_from_concentration(8.0).unwrap();
        assert!((alpha - 8.0 / 4.1).abs() < 1e-10);
    }

    #[test]
    fn test_mcbride_z0_equals_m0() {
        let m = mass_accretion_mcbride(1e12, 0.0, 0.1, 0.8).unwrap();
        assert!((m - 1e12).abs() < 1.0);
    }

    #[test]
    fn test_mcbride_decreases_with_z() {
        let m0 = mass_accretion_mcbride(1e12, 0.0, 0.1, 0.8).unwrap();
        let m1 = mass_accretion_mcbride(1e12, 1.0, 0.1, 0.8).unwrap();
        let m3 = mass_accretion_mcbride(1e12, 3.0, 0.1, 0.8).unwrap();
        assert!(m0 > m1 && m1 > m3);
    }

    #[test]
    fn test_mah_invalid_inputs() {
        assert!(mass_accretion_wechsler(0.0, 0.0, 1.0).is_err());
        assert!(mass_accretion_wechsler(1e12, 0.0, 0.0).is_err());
        assert!(mass_accretion_wechsler(f64::NAN, 0.0, 1.0).is_err());
        assert!(mass_accretion_mcbride(0.0, 0.0, 0.1, 0.8).is_err());
        assert!(mass_accretion_mcbride(1e12, 0.0, 0.1, 0.0).is_err());
    }

    // -- SHAM / Schechter --

    #[test]
    fn test_schechter_smf_positive() {
        let phi = schechter_smf(10.0, 10.65, 4.5e-3, -1.2).unwrap();
        assert!(phi > 0.0);
    }

    #[test]
    fn test_schechter_smf_cutoff() {
        let phi_below = schechter_smf(10.0, 10.65, 4.5e-3, -1.2).unwrap();
        let phi_above = schechter_smf(12.0, 10.65, 4.5e-3, -1.2).unwrap();
        assert!(phi_below > phi_above, "exponential cutoff above M*₀");
    }

    #[test]
    fn test_schechter_cumulative_monotonic() {
        let n1 = schechter_cumulative(9.0, 10.65, 4.5e-3, -1.2).unwrap();
        let n2 = schechter_cumulative(10.0, 10.65, 4.5e-3, -1.2).unwrap();
        let n3 = schechter_cumulative(11.0, 10.65, 4.5e-3, -1.2).unwrap();
        assert!(n1 > n2 && n2 > n3, "cumulative should decrease");
    }

    #[test]
    fn test_sham_stellar_mass_reasonable() {
        let c = cosmo();
        let log_ms = sham_stellar_mass(1e12, 0.0, 10.65, 4.5e-3, -1.2, &c).unwrap();
        assert!(
            log_ms > 9.0 && log_ms < 12.0,
            "MW halo stellar mass log M* = {log_ms}"
        );
    }

    #[test]
    fn test_sham_increases_with_halo_mass() {
        let c = cosmo();
        let ms1 = sham_stellar_mass(1e11, 0.0, 10.65, 4.5e-3, -1.2, &c).unwrap();
        let ms2 = sham_stellar_mass(1e12, 0.0, 10.65, 4.5e-3, -1.2, &c).unwrap();
        let ms3 = sham_stellar_mass(1e13, 0.0, 10.65, 4.5e-3, -1.2, &c).unwrap();
        assert!(ms3 > ms2 && ms2 > ms1, "more massive halos → more stars");
    }

    #[test]
    fn test_sham_invalid() {
        let c = cosmo();
        assert!(sham_stellar_mass(0.0, 0.0, 10.65, 4.5e-3, -1.2, &c).is_err());
        assert!(sham_stellar_mass(f64::NAN, 0.0, 10.65, 4.5e-3, -1.2, &c).is_err());
    }
}
