//! Matter power spectrum — structure formation, transfer functions.
//!
//! The power spectrum P(k) describes the amplitude of density fluctuations
//! as a function of wavenumber k. Combined with the growth factor, it
//! predicts structure at all scales and redshifts.
//!
//! All functions take a [`Cosmology`] reference instead of individual
//! cosmological parameters.

use crate::constants::C;
use crate::cosmology::{Cosmology, FilterFunction};
use crate::error::{BrahmandaError, ensure_finite, require_finite};

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
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::transfer_function_eh;
///
/// let cosmo = Cosmology::planck2018();
/// let t = transfer_function_eh(0.001, cosmo.omega_m, cosmo.omega_b, cosmo.h).unwrap();
/// assert!(t > 0.8);
/// ```
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
    let s = 44.5
        * (omega_mh2).sqrt().recip()
        * (1.0 + (omega_bh2 / 0.024).powf(0.46)).recip()
        * (omega_mh2).powf(0.25)
        * 100.0; // rough approximation in Mpc/h

    // Silk damping scale
    let alpha_gamma = 1.0 - 0.328 * (omega_bh2 / omega_mh2).ln() * (omega_b / omega_m)
        + 0.38 * (omega_bh2 / omega_mh2).ln() * (omega_b / omega_m).powi(2);

    // Shape parameter
    let gamma_eff =
        omega_m * h * (alpha_gamma + (1.0 - alpha_gamma) / (1.0 + (0.43 * k_mpc * s).powi(4)));

    // BBKS transfer function with effective shape
    let q = k_mpc * theta * theta / gamma_eff;
    let t = (1.0 + 3.89 * q + (16.1 * q).powi(2) + (5.46 * q).powi(3) + (6.71 * q).powi(4))
        .powf(-0.25)
        * (2.34 * q).ln_1p()
        / (2.34 * q);

    ensure_finite(t, "transfer_function_eh")
}

/// Eisenstein-Hu transfer function T(k) — full version with BAO wiggles.
///
/// Implements the baryon + CDM transfer function from Eisenstein & Hu (1998)
/// including baryon acoustic oscillations.
///
/// # Arguments
/// * `k_mpc` — Wavenumber in h/Mpc.
/// * `omega_m` — Matter density parameter.
/// * `omega_b` — Baryon density parameter.
/// * `h` — Dimensionless Hubble parameter (H₀/100).
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::transfer_function_eh_wiggle;
///
/// let cosmo = Cosmology::planck2018();
/// let t = transfer_function_eh_wiggle(0.001, cosmo.omega_m, cosmo.omega_b, cosmo.h).unwrap();
/// assert!(t > 0.8);
/// ```
pub fn transfer_function_eh_wiggle(
    k_mpc: f64,
    omega_m: f64,
    omega_b: f64,
    h: f64,
) -> Result<f64, BrahmandaError> {
    require_finite(k_mpc, "transfer_function_eh_wiggle")?;
    require_finite(omega_m, "transfer_function_eh_wiggle")?;
    require_finite(omega_b, "transfer_function_eh_wiggle")?;
    require_finite(h, "transfer_function_eh_wiggle")?;
    if k_mpc <= 0.0 || omega_m <= 0.0 || h <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "transfer_function_eh_wiggle: k, omega_m, h must be positive".to_string(),
        ));
    }
    if omega_b < 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "transfer_function_eh_wiggle: omega_b must be non-negative".to_string(),
        ));
    }

    let omega_mh2 = omega_m * h * h;
    let omega_bh2 = omega_b * h * h;
    let f_b = omega_b / omega_m;
    let f_c = 1.0 - f_b;
    let theta_27: f64 = 2.725 / 2.7; // T_CMB / 2.7 K

    // Eq. 2: redshift of matter-radiation equality
    let z_eq = 2.5e4 * omega_mh2 * theta_27.powi(-4);
    let k_eq = 7.46e-2 * omega_mh2 * theta_27.powi(-2); // Eq. 3: in h/Mpc

    // Eq. 4: redshift of baryon drag epoch
    let b1_d = 0.313 * omega_mh2.powf(-0.419) * (1.0 + 0.607 * omega_mh2.powf(0.674));
    let b2_d = 0.238 * omega_mh2.powf(0.223);
    let z_d = 1291.0 * omega_mh2.powf(0.251) / (1.0 + 0.659 * omega_mh2.powf(0.828))
        * (1.0 + b1_d * omega_bh2.powf(b2_d));

    // Eq. 5: baryon-to-photon ratio R = 3ρ_b/(4ρ_γ) at z_eq and z_d
    let r_coeff = 31.5e3 * omega_bh2 * theta_27.powi(-4); // in units of 1/z
    let r_eq = r_coeff / z_eq;
    let r_d = r_coeff / z_d;

    // Eq. 6: sound horizon at drag epoch (Mpc/h)
    let s = (2.0 / (3.0 * k_eq))
        * (6.0 / r_eq).sqrt()
        * (((1.0 + r_d).sqrt() + (r_d + r_eq).sqrt()) / (1.0 + r_eq.sqrt())).ln();

    // Eq. 7: Silk damping scale (1/Mpc h)
    let k_silk =
        1.6 * omega_bh2.powf(0.52) * omega_mh2.powf(0.73) * (1.0 + (10.4 * omega_mh2).powf(-0.95));

    // Eq. 10: CDM alpha_c
    let a1 = (46.9 * omega_mh2).powf(0.670) * (1.0 + (32.1 * omega_mh2).powf(-0.532));
    let a2 = (12.0 * omega_mh2).powf(0.424) * (1.0 + (45.0 * omega_mh2).powf(-0.582));
    let alpha_c = a1.powf(-f_b) * a2.powf(-f_b.powi(3));

    // Eq. 12: CDM beta_c
    let b1_c = 0.944 / (1.0 + (458.0 * omega_mh2).powf(-0.708));
    let b2_c = (0.395 * omega_mh2).powf(-0.0266);
    let beta_c = 1.0 / (1.0 + b1_c * (f_c.powf(b2_c) - 1.0));

    // Eq. 10 (normalized wavenumber)
    let q = k_mpc / (13.41 * k_eq);

    // Eq. 17-18: T_0(k|q, αc, βc) — CDM transfer function form
    let t0 = |q_val: f64, ac: f64, bc: f64| -> f64 {
        let c_val = 14.2 / ac + 386.0 / (1.0 + 69.9 * q_val.powf(1.08));
        let l = (std::f64::consts::E + 1.8 * bc * q_val).ln();
        l / (l + c_val * q_val * q_val)
    };

    // Eq. 17: CDM piece with interpolation
    let f_val = 1.0 / (1.0 + (k_mpc * s / 5.4).powi(4));
    let t_c = f_val * t0(q, 1.0, beta_c) + (1.0 - f_val) * t0(q, alpha_c, 1.0);

    // Eq. 8: G(y) function for baryon suppression
    let y = z_eq / z_d;
    let sqrt_1py = (1.0 + y).sqrt();
    let g_y = y * (-6.0 * sqrt_1py + (2.0 + 3.0 * y) * ((sqrt_1py + 1.0) / (sqrt_1py - 1.0)).ln());

    // Eq. 14: baryon alpha_b
    let alpha_b = 2.07 * k_eq * s * (1.0 + r_d).powf(-3.0 / 4.0) * g_y;

    // Eq. 24: baryon beta_node
    let beta_node = 8.41 * omega_mh2.powf(0.435);

    // Eq. 23: baryon beta_b
    let beta_b = 0.5 + f_b + (3.0 - 2.0 * f_b) * ((17.2 * omega_mh2).powi(2) + 1.0).sqrt();

    // Eq. 22: effective sound horizon with node shift
    let ks = k_mpc * s;
    let s_tilde = s / (1.0 + (beta_node / ks).powi(3)).cbrt();

    // Eq. 19-21: baryon transfer function
    let ks_tilde = k_mpc * s_tilde;
    let j0_tilde = if ks_tilde.abs() < 1e-10 {
        1.0
    } else {
        ks_tilde.sin() / ks_tilde
    };

    let t0_b_val = t0(q, 1.0, 1.0);
    let t_b = (t0_b_val / (1.0 + (ks / 5.2).powi(2))
        + alpha_b / (1.0 + (beta_b / ks).powi(3)) * (-(k_mpc / k_silk).powf(1.4)).exp())
        * j0_tilde;

    // Eq. 16: full transfer function
    let t = f_b * t_b + f_c * t_c;
    ensure_finite(t.abs(), "transfer_function_eh_wiggle")
}

/// Primordial power spectrum: P_prim(k) ∝ k^n_s.
///
/// Returns the unnormalized primordial power at wavenumber k.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::primordial_power;
///
/// let cosmo = Cosmology::planck2018();
/// let p = primordial_power(1.0, &cosmo).unwrap();
/// assert!((p - 1.0).abs() < 1e-10);
/// ```
#[inline]
pub fn primordial_power(k_mpc: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(k_mpc, "primordial_power")?;
    if k_mpc <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "primordial_power: k must be positive".to_string(),
        ));
    }
    ensure_finite(k_mpc.powf(cosmo.n_s), "primordial_power")
}

/// Linear growth factor D(z) — Carroll, Press & Turner (1992) approximation.
///
/// D(z) / D(0) for a flat ΛCDM cosmology.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::growth_factor;
///
/// let cosmo = Cosmology::planck2018();
/// let d = growth_factor(0.0, &cosmo).unwrap();
/// assert!((d - 1.0).abs() < 1e-10);
///
/// let d1 = growth_factor(1.0, &cosmo).unwrap();
/// assert!(d1 < d);
/// ```
pub fn growth_factor(z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(z, "growth_factor")?;
    if z < -1.0 {
        return Err(BrahmandaError::InvalidStructure(
            "growth_factor: z must be >= -1".to_string(),
        ));
    }

    let omega_m = cosmo.omega_m;
    let omega_lambda = cosmo.omega_lambda;

    let growth_at = |zz: f64| -> f64 {
        let a = 1.0 / (1.0 + zz);
        let om_z = omega_m / (omega_m + omega_lambda * a * a * a);
        let ol_z = omega_lambda * a * a * a / (omega_m + omega_lambda * a * a * a);
        // CPT92 approximation
        let d =
            2.5 * om_z / (om_z.powf(4.0 / 7.0) - ol_z + (1.0 + om_z / 2.0) * (1.0 + ol_z / 70.0));
        d * a
    };

    let d_z = growth_at(z);
    let d_0 = growth_at(0.0);
    ensure_finite(d_z / d_0, "growth_factor")
}

/// Linear power spectrum P_lin(k, z).
///
/// P_lin(k, z) = A_s × k^n_s × T²(k) × D²(z)
///
/// where A_s is the amplitude normalized to σ_8.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::linear_power_spectrum;
///
/// let cosmo = Cosmology::planck2018();
/// let p = linear_power_spectrum(0.1, 0.0, &cosmo).unwrap();
/// assert!(p > 0.0);
/// ```
pub fn linear_power_spectrum(k_mpc: f64, z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(k_mpc, "linear_power_spectrum")?;
    require_finite(z, "linear_power_spectrum")?;
    if k_mpc <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "linear_power_spectrum: k must be positive".to_string(),
        ));
    }
    let a_s = cosmo.power_spectrum_amplitude();
    let t = transfer_function_eh(k_mpc, cosmo.omega_m, cosmo.omega_b, cosmo.h)?;
    let d = growth_factor(z, cosmo)?;
    ensure_finite(
        a_s * k_mpc.powf(cosmo.n_s) * t * t * d * d,
        "linear_power_spectrum",
    )
}

/// σ(R, z) — RMS mass fluctuation in a top-hat sphere of radius R (Mpc/h).
///
/// Exact computation via adaptive Simpson integration of the power spectrum
/// with a real-space top-hat window function.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::sigma_r;
///
/// let cosmo = Cosmology::planck2018();
/// let s = sigma_r(8.0, 0.0, &cosmo).unwrap();
/// assert!((s - 0.811).abs() < 0.01, "σ(8,0) = {s}");
/// ```
pub fn sigma_r(r_mpc: f64, z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    sigma_r_filter(r_mpc, z, cosmo, FilterFunction::TopHat)
}

/// σ(R, z) with a specified filter function.
///
/// Computes σ²(R) = A_s ∫ k³ k^n_s T²(k) W²(kR) dk / (2π²) in log-k space,
/// then multiplies by D²(z) and takes the square root.
///
/// ```
/// use brahmanda::{Cosmology, FilterFunction};
/// use brahmanda::power_spectrum::sigma_r_filter;
///
/// let cosmo = Cosmology::planck2018();
/// let s = sigma_r_filter(8.0, 0.0, &cosmo, FilterFunction::TopHat).unwrap();
/// assert!((s - cosmo.sigma_8).abs() < 0.01);
/// ```
pub fn sigma_r_filter(
    r_mpc: f64,
    z: f64,
    cosmo: &Cosmology,
    filter: FilterFunction,
) -> Result<f64, BrahmandaError> {
    require_finite(r_mpc, "sigma_r_filter")?;
    require_finite(z, "sigma_r_filter")?;
    if r_mpc <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "sigma_r_filter: radius must be positive".to_string(),
        ));
    }

    let a_s = cosmo.power_spectrum_amplitude();
    let n_s = cosmo.n_s;
    let omega_m = cosmo.omega_m;
    let omega_b = cosmo.omega_b;
    let h = cosmo.h;

    let ln_k_min = -4.0_f64 * std::f64::consts::LN_10; // ln(1e-4)
    let ln_k_max = (100.0_f64).ln();

    let integrand = |u: f64| -> f64 {
        let k = u.exp();
        let t = transfer_function_eh(k, omega_m, omega_b, h).unwrap_or(0.0);
        let w = filter.window(k * r_mpc);
        k.powi(3) * k.powf(n_s) * t * t * w * w
            / (2.0 * std::f64::consts::PI * std::f64::consts::PI)
    };

    let sigma_sq_raw = hisab::calc::integral_adaptive_simpson(integrand, ln_k_min, ln_k_max, 1e-8)
        .map_err(|e| BrahmandaError::Computation(format!("sigma_r integration: {e}")))?;

    let d = growth_factor(z, cosmo)?;
    let sigma = (a_s * sigma_sq_raw).sqrt() * d;
    ensure_finite(sigma, "sigma_r_filter")
}

/// Dark energy equation of state — CPL parameterization.
///
/// w(z) = w₀ + w_a × z/(1+z)
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::dark_energy_eos;
///
/// let cosmo = Cosmology::planck2018();
/// let w = dark_energy_eos(0.5, &cosmo).unwrap();
/// assert!((w - (-1.0)).abs() < 1e-10);
/// ```
#[inline]
pub fn dark_energy_eos(z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(z, "dark_energy_eos")?;
    if z < -1.0 {
        return Err(BrahmandaError::InvalidStructure(
            "dark_energy_eos: z must be >= -1".to_string(),
        ));
    }
    ensure_finite(cosmo.w0 + cosmo.wa * z / (1.0 + z), "dark_energy_eos")
}

/// Dark energy density evolution factor for CPL parameterization.
///
/// ρ_DE(z)/ρ_DE(0) = (1+z)^(3(1+w₀+w_a)) × exp(-3 w_a z/(1+z))
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::dark_energy_density_ratio;
///
/// let cosmo = Cosmology::planck2018();
/// let ratio = dark_energy_density_ratio(1.0, &cosmo).unwrap();
/// assert!((ratio - 1.0).abs() < 1e-10); // ΛCDM
/// ```
pub fn dark_energy_density_ratio(z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(z, "dark_energy_density_ratio")?;
    if z < -1.0 {
        return Err(BrahmandaError::InvalidStructure(
            "dark_energy_density_ratio: z must be >= -1".to_string(),
        ));
    }
    let a = 1.0 / (1.0 + z);
    let exponent = 3.0 * (1.0 + cosmo.w0 + cosmo.wa);
    let ratio = (1.0 + z).powf(exponent) * (-3.0 * cosmo.wa * (1.0 - a)).exp();
    ensure_finite(ratio, "dark_energy_density_ratio")
}

/// Hubble parameter E(z) = H(z)/H₀ for w₀w_a dark energy.
///
/// E²(z) = Ω_r(1+z)⁴ + Ω_m(1+z)³ + Ω_k(1+z)² + Ω_DE × (ρ_DE(z)/ρ_DE(0))
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::hubble_parameter_ratio;
///
/// let cosmo = Cosmology::planck2018();
/// let e = hubble_parameter_ratio(0.0, &cosmo).unwrap();
/// assert!((e - 1.0).abs() < 1e-6);
/// ```
pub fn hubble_parameter_ratio(z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(z, "hubble_parameter_ratio")?;
    if z < -1.0 {
        return Err(BrahmandaError::InvalidStructure(
            "hubble_parameter_ratio: z must be >= -1".to_string(),
        ));
    }
    let opz = 1.0 + z;
    let de_ratio = dark_energy_density_ratio(z, cosmo)?;
    let e2 = cosmo.omega_r * opz.powi(4)
        + cosmo.omega_m * opz.powi(3)
        + cosmo.omega_k * opz.powi(2)
        + cosmo.omega_lambda * de_ratio;
    ensure_finite(e2.sqrt(), "hubble_parameter_ratio")
}

/// Linear growth factor for w₀w_a dark energy — CPT92 fitting formula.
///
/// Returns D(z)/D(0).
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::growth_factor_w;
///
/// let cosmo = Cosmology::planck2018();
/// let d = growth_factor_w(1.0, &cosmo).unwrap();
/// assert!(d > 0.0 && d < 1.0);
/// ```
pub fn growth_factor_w(z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(z, "growth_factor_w")?;
    if z < -1.0 {
        return Err(BrahmandaError::InvalidStructure(
            "growth_factor_w: z must be >= -1".to_string(),
        ));
    }

    let growth_at = |zz: f64| -> Result<f64, BrahmandaError> {
        let a = 1.0 / (1.0 + zz);
        let e = hubble_parameter_ratio(zz, cosmo)?;
        let e2 = e * e;
        let om_z = cosmo.omega_m * (1.0 + zz).powi(3) / e2;
        let ol_z = 1.0 - om_z;
        let d =
            2.5 * om_z / (om_z.powf(4.0 / 7.0) - ol_z + (1.0 + om_z / 2.0) * (1.0 + ol_z / 70.0));
        Ok(d * a)
    };

    let d_z = growth_at(z)?;
    let d_0 = growth_at(0.0)?;
    ensure_finite(d_z / d_0, "growth_factor_w")
}

/// Comoving distance χ(z) in Mpc (physical, not h⁻¹ Mpc).
///
/// χ(z) = (c/H₀) ∫₀ᶻ dz'/E(z')
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::comoving_distance;
///
/// let cosmo = Cosmology::planck2018();
/// let chi = comoving_distance(0.0, &cosmo).unwrap();
/// assert!(chi.abs() < 1e-10);
///
/// let chi = comoving_distance(1.0, &cosmo).unwrap();
/// assert!(chi > 2000.0 && chi < 5000.0);
/// ```
pub fn comoving_distance(z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(z, "comoving_distance")?;
    if z < 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "comoving_distance: z must be non-negative".to_string(),
        ));
    }
    if z == 0.0 {
        return Ok(0.0);
    }

    let n = 200_usize;
    let dz = z / n as f64;
    let mut sum = 0.0;

    for i in 0..n {
        let z0 = i as f64 * dz;
        let z1 = z0 + dz / 2.0;
        let z2 = z0 + dz;
        let e0 = hubble_parameter_ratio(z0, cosmo)?;
        let e1 = hubble_parameter_ratio(z1, cosmo)?;
        let e2 = hubble_parameter_ratio(z2, cosmo)?;
        sum += (dz / 6.0) * (1.0 / e0 + 4.0 / e1 + 1.0 / e2);
    }

    let h0_km = cosmo.h * 100.0;
    let c_over_h0 = (C / 1e3) / h0_km;
    ensure_finite(c_over_h0 * sum, "comoving_distance")
}

/// Angular power spectrum C_l — Limber approximation.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::angular_power_spectrum_limber;
///
/// let cosmo = Cosmology::planck2018();
/// let cl = angular_power_spectrum_limber(10, 1.0, &cosmo).unwrap();
/// assert!(cl > 0.0);
/// ```
pub fn angular_power_spectrum_limber(
    l: u32,
    z_max: f64,
    cosmo: &Cosmology,
) -> Result<f64, BrahmandaError> {
    require_finite(z_max, "angular_power_spectrum_limber")?;
    if l == 0 {
        return Err(BrahmandaError::InvalidStructure(
            "angular_power_spectrum_limber: l must be > 0".to_string(),
        ));
    }
    if z_max <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "angular_power_spectrum_limber: z_max must be positive".to_string(),
        ));
    }

    let l_f = l as f64;
    let chi_max = comoving_distance(z_max, cosmo)?;

    let n = 100_usize;
    let dz = z_max / n as f64;
    let mut sum = 0.0;

    let integrand = |z: f64| -> Result<f64, BrahmandaError> {
        if z < 1e-6 {
            return Ok(0.0);
        }
        let chi = comoving_distance(z, cosmo)?;
        if chi < 1e-6 {
            return Ok(0.0);
        }
        let k = (l_f + 0.5) / chi;
        let k_h = k / cosmo.h;

        if k_h <= 0.0 || k_h > 100.0 {
            return Ok(0.0);
        }

        let t = transfer_function_eh(k_h, cosmo.omega_m, cosmo.omega_b, cosmo.h)?;
        let d = growth_factor_w(z, cosmo)?;

        let pk = k_h.powf(cosmo.n_s) * t * t * d * d;

        let e = hubble_parameter_ratio(z, cosmo)?;

        Ok(pk / (chi * chi) / e)
    };

    for i in 0..n {
        let z0 = i as f64 * dz;
        let z1 = z0 + dz / 2.0;
        let z2 = z0 + dz;
        let f0 = integrand(z0)?;
        let f1 = integrand(z1)?;
        let f2 = integrand(z2)?;
        sum += (dz / 6.0) * (f0 + 4.0 * f1 + f2);
    }

    let h0_km = cosmo.h * 100.0;
    let c_over_h0 = (C / 1e3) / h0_km;
    let cl = c_over_h0 * sum / (chi_max * chi_max);
    ensure_finite(cl, "angular_power_spectrum_limber")
}

/// Luminosity distance d_L(z) in Mpc.
///
/// d_L = (1+z) × χ(z) for a flat universe.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::luminosity_distance;
///
/// let cosmo = Cosmology::planck2018();
/// let dl = luminosity_distance(0.0, &cosmo).unwrap();
/// assert!(dl.abs() < 1e-10);
/// ```
pub fn luminosity_distance(z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    let chi = comoving_distance(z, cosmo)?;
    ensure_finite((1.0 + z) * chi, "luminosity_distance")
}

/// Distance modulus μ(z) for SN Ia cosmology.
///
/// μ = 5 log₁₀(d_L \[Mpc\]) + 25
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::distance_modulus;
///
/// let cosmo = Cosmology::planck2018();
/// let mu = distance_modulus(1.0, &cosmo).unwrap();
/// assert!(mu > 42.0 && mu < 46.0);
/// ```
pub fn distance_modulus(z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(z, "distance_modulus")?;
    if z <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "distance_modulus: z must be positive".to_string(),
        ));
    }
    let dl = luminosity_distance(z, cosmo)?;
    ensure_finite(5.0 * dl.log10() + 25.0, "distance_modulus")
}

/// Angular diameter distance d_A(z) in Mpc.
///
/// d_A = χ(z) / (1+z) for a flat universe.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::angular_diameter_distance;
///
/// let cosmo = Cosmology::planck2018();
/// let da = angular_diameter_distance(1.0, &cosmo).unwrap();
/// assert!(da > 1000.0 && da < 2500.0);
/// ```
pub fn angular_diameter_distance(z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    let chi = comoving_distance(z, cosmo)?;
    ensure_finite(chi / (1.0 + z), "angular_diameter_distance")
}

/// Comoving volume element dV/dz per steradian (Mpc³/sr).
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::comoving_volume_element;
///
/// let cosmo = Cosmology::planck2018();
/// let dv = comoving_volume_element(1.0, &cosmo).unwrap();
/// assert!(dv > 0.0);
/// ```
pub fn comoving_volume_element(z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(z, "comoving_volume_element")?;
    if z < 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "comoving_volume_element: z must be non-negative".to_string(),
        ));
    }
    let chi = comoving_distance(z, cosmo)?;
    let e = hubble_parameter_ratio(z, cosmo)?;
    let h0_km = cosmo.h * 100.0;
    let c_over_h0 = (C / 1e3) / h0_km;
    ensure_finite(c_over_h0 * chi * chi / e, "comoving_volume_element")
}

/// Lookback time t_lb(z) in Gyr.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::lookback_time;
///
/// let cosmo = Cosmology::planck2018();
/// let t = lookback_time(1.0, &cosmo).unwrap();
/// assert!(t > 6.0 && t < 10.0);
/// ```
pub fn lookback_time(z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(z, "lookback_time")?;
    if z < 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "lookback_time: z must be non-negative".to_string(),
        ));
    }
    if z == 0.0 {
        return Ok(0.0);
    }

    let n = 200_usize;
    let ln1z = (1.0 + z).ln();
    let d_ln = ln1z / n as f64;
    let mut sum = 0.0;

    for i in 0..n {
        let u0 = i as f64 * d_ln;
        let u1 = u0 + d_ln / 2.0;
        let u2 = u0 + d_ln;
        let z0 = u0.exp() - 1.0;
        let z1 = u1.exp() - 1.0;
        let z2 = u2.exp() - 1.0;
        let e0 = hubble_parameter_ratio(z0, cosmo)?;
        let e1 = hubble_parameter_ratio(z1, cosmo)?;
        let e2 = hubble_parameter_ratio(z2, cosmo)?;
        sum += (d_ln / 6.0) * (1.0 / e0 + 4.0 / e1 + 1.0 / e2);
    }

    let h0_inv_gyr = 1.0 / (cosmo.h0_si * 3.15576e16);
    ensure_finite(h0_inv_gyr * sum, "lookback_time")
}

/// Age of the universe at redshift z in Gyr.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::age_of_universe;
///
/// let cosmo = Cosmology::planck2018();
/// let age = age_of_universe(0.0, &cosmo).unwrap();
/// assert!(age > 12.0 && age < 15.0);
/// ```
pub fn age_of_universe(z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(z, "age_of_universe")?;
    if z < 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "age_of_universe: z must be non-negative".to_string(),
        ));
    }

    let z_max = 2000.0;
    if z >= z_max {
        return Ok(0.0);
    }

    let n = 600_usize;
    let ln1z_lo = (1.0 + z).ln();
    let ln1z_hi = (1.0 + z_max).ln();
    let d_ln = (ln1z_hi - ln1z_lo) / n as f64;
    let mut sum = 0.0;

    for i in 0..n {
        let u0 = ln1z_lo + i as f64 * d_ln;
        let u1 = u0 + d_ln / 2.0;
        let u2 = u0 + d_ln;
        let z0 = u0.exp() - 1.0;
        let z1 = u1.exp() - 1.0;
        let z2 = u2.exp() - 1.0;
        let e0 = hubble_parameter_ratio(z0, cosmo)?;
        let e1 = hubble_parameter_ratio(z1, cosmo)?;
        let e2 = hubble_parameter_ratio(z2, cosmo)?;
        sum += (d_ln / 6.0) * (1.0 / e0 + 4.0 / e1 + 1.0 / e2);
    }

    let h0_inv_gyr = 1.0 / (cosmo.h0_si * 3.15576e16);
    ensure_finite(h0_inv_gyr * sum, "age_of_universe")
}

/// Linear growth rate f(z) = d ln D / d ln a.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::growth_rate;
///
/// let cosmo = Cosmology::planck2018();
/// let f = growth_rate(0.0, &cosmo).unwrap();
/// assert!(f > 0.4 && f < 0.7);
/// ```
pub fn growth_rate(z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(z, "growth_rate")?;
    if z < -1.0 {
        return Err(BrahmandaError::InvalidStructure(
            "growth_rate: z must be >= -1".to_string(),
        ));
    }

    let eps = 1e-4;
    let d_hi = growth_factor_w(z + eps, cosmo)?;
    let d_lo = growth_factor_w((z - eps).max(0.0), cosmo)?;
    let dz_eff = if z < eps { z + eps } else { 2.0 * eps };
    let dlnd_dz = (d_hi.ln() - d_lo.ln()) / dz_eff;
    let f = -(1.0 + z) * dlnd_dz;
    ensure_finite(f, "growth_rate")
}

/// Effective spectral index n_eff at scale R (Mpc/h).
///
/// n_eff = -2 × d ln σ(R) / d ln R - 3
fn effective_spectral_index(r_mpc: f64, z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    let dr = r_mpc * 0.01;
    let s_lo = sigma_r(r_mpc - dr, z, cosmo)?;
    let s_hi = sigma_r(r_mpc + dr, z, cosmo)?;
    let dlns_dlnr = (s_hi.ln() - s_lo.ln()) / ((r_mpc + dr).ln() - (r_mpc - dr).ln());
    ensure_finite(-2.0 * dlns_dlnr - 3.0, "effective_spectral_index")
}

/// Nonlinear power spectrum — Smith et al. (2003) Halofit.
///
/// Returns the nonlinear dimensionless power spectrum Δ²_NL(k) = k³ P_NL(k) / (2π²).
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::halofit_nl;
///
/// let cosmo = Cosmology::planck2018();
/// let dnl = halofit_nl(0.1, 0.0, &cosmo).unwrap();
/// assert!(dnl > 0.0);
/// ```
pub fn halofit_nl(k_mpc: f64, z: f64, cosmo: &Cosmology) -> Result<f64, BrahmandaError> {
    require_finite(k_mpc, "halofit_nl")?;
    require_finite(z, "halofit_nl")?;
    if k_mpc <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "halofit_nl: k must be positive".to_string(),
        ));
    }

    // Find nonlinear scale: R_nl where σ(R_nl, z) = 1
    let mut r_lo = 0.01_f64;
    let mut r_hi = 100.0_f64;
    for _ in 0..60 {
        let r_mid = (r_lo * r_hi).sqrt();
        let s = sigma_r(r_mid, z, cosmo)?;
        if s > 1.0 {
            r_lo = r_mid;
        } else {
            r_hi = r_mid;
        }
    }
    let r_nl = (r_lo * r_hi).sqrt();
    let k_sigma = 1.0 / r_nl;

    let n_eff = effective_spectral_index(r_nl, z, cosmo)?;
    let n2 = n_eff * n_eff;

    let y = k_mpc / k_sigma;

    // Smith et al. 2003, Table 2
    let a_n = 10.0_f64.powf(
        1.4861 + 1.8369 * n_eff + 1.6762 * n2 + 0.7940 * n_eff * n2 + 0.1670 * n2 * n2
            - 0.6206 * n_eff.powi(5),
    );
    let b_n = 10.0_f64.powf(0.9463 + 0.9466 * n_eff + 0.3084 * n2 - 0.9400 * n_eff * n2);
    let c_n = 10.0_f64.powf(-0.2807 + 0.6669 * n_eff + 0.3214 * n2 - 0.0793 * n_eff * n2);
    let mu_n = 10.0_f64.powf(-3.5442 + 0.1908 * n_eff);
    let nu_n = 10.0_f64.powf(0.9589 + 1.2857 * n_eff);
    let alpha_n = (1.3884 + 0.3700 * n_eff - 0.1452 * n2).abs();
    let beta_n = 0.8291 + 0.9854 * n_eff + 0.3401 * n2;

    // Linear Δ²(k) ∝ σ²(1/k) approximately
    let sigma_k = sigma_r(1.0 / k_mpc, z, cosmo)?;
    let delta2_lin = sigma_k * sigma_k;

    // Two-halo (quasi-linear) term: Δ²_Q
    let f_y = y / 4.0 + y * y / 8.0;
    let delta2_q = delta2_lin
        * ((1.0 + delta2_lin).powf(beta_n) / (1.0 + alpha_n * delta2_lin))
        * (-f_y).exp();

    // One-halo term: Δ²_H
    let delta2_h =
        a_n * y.powf(3.0 * mu_n) / (1.0 + b_n * y.powf(-1.0)) / (1.0 + c_n * y.powf(3.0 * nu_n));

    ensure_finite(delta2_q + delta2_h, "halofit_nl")
}

/// Ordinary Sachs-Wolfe effect: ΔT/T from gravitational potential.
///
/// (ΔT/T)_SW = Φ / 3
///
/// ```
/// use brahmanda::power_spectrum::sachs_wolfe;
///
/// let dt_t = sachs_wolfe(1e-5).unwrap();
/// assert!((dt_t - 1e-5 / 3.0).abs() < 1e-15);
/// ```
#[inline]
pub fn sachs_wolfe(phi: f64) -> Result<f64, BrahmandaError> {
    require_finite(phi, "sachs_wolfe")?;
    ensure_finite(phi / 3.0, "sachs_wolfe")
}

/// Integrated Sachs-Wolfe effect: ΔT/T from time-varying potentials.
///
/// ```
/// use brahmanda::Cosmology;
/// use brahmanda::power_spectrum::integrated_sachs_wolfe;
///
/// let cosmo = Cosmology::planck2018();
/// let isw = integrated_sachs_wolfe(0.0, 2.0, &cosmo).unwrap();
/// assert!(isw.abs() > 0.0);
/// ```
pub fn integrated_sachs_wolfe(
    z_min: f64,
    z_max: f64,
    cosmo: &Cosmology,
) -> Result<f64, BrahmandaError> {
    require_finite(z_min, "integrated_sachs_wolfe")?;
    require_finite(z_max, "integrated_sachs_wolfe")?;
    if z_min < 0.0 || z_max <= z_min {
        return Err(BrahmandaError::InvalidStructure(
            "integrated_sachs_wolfe: need 0 <= z_min < z_max".to_string(),
        ));
    }

    let n = 200_usize;
    let dz = (z_max - z_min) / n as f64;
    let eps = 1e-4;
    let mut sum = 0.0;

    let isw_integrand = |z: f64| -> Result<f64, BrahmandaError> {
        let d = growth_factor_w(z, cosmo)?;
        let e = hubble_parameter_ratio(z, cosmo)?;

        let d_hi = growth_factor_w(z + eps, cosmo)?;
        let d_lo = growth_factor_w((z - eps).max(0.0), cosmo)?;
        let dz_eff = if z < eps { z + eps } else { 2.0 * eps };
        let dlnd_dz = (d_hi.ln() - d_lo.ln()) / dz_eff;
        let f_growth = -(1.0 + z) * dlnd_dz;

        Ok((f_growth - 1.0) * d * e)
    };

    for i in 0..n {
        let z0 = z_min + i as f64 * dz;
        let z1 = z0 + dz / 2.0;
        let z2 = z0 + dz;
        let f0 = isw_integrand(z0)?;
        let f1 = isw_integrand(z1)?;
        let f2 = isw_integrand(z2)?;
        sum += (dz / 6.0) * (f0 + 4.0 * f1 + f2);
    }

    ensure_finite(2.0 * sum, "integrated_sachs_wolfe")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Cosmology;

    fn cosmo() -> Cosmology {
        Cosmology::planck2018()
    }

    #[test]
    fn test_transfer_function_small_k() {
        let c = cosmo();
        let t = transfer_function_eh(0.001, c.omega_m, c.omega_b, c.h).unwrap();
        assert!(t > 0.8, "T(k→0) should approach 1, got {t}");
    }

    #[test]
    fn test_transfer_function_large_k() {
        let c = cosmo();
        let t = transfer_function_eh(1.0, c.omega_m, c.omega_b, c.h).unwrap();
        assert!(t < 0.5, "T(k=1) should be suppressed, got {t}");
    }

    #[test]
    fn test_growth_factor_z0() {
        let c = cosmo();
        let d = growth_factor(0.0, &c).unwrap();
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_growth_factor_decreases_with_z() {
        let c = cosmo();
        let d0 = growth_factor(0.0, &c).unwrap();
        let d1 = growth_factor(1.0, &c).unwrap();
        let d5 = growth_factor(5.0, &c).unwrap();
        assert!(d0 > d1);
        assert!(d1 > d5);
    }

    #[test]
    fn test_sigma_8_at_z0() {
        let c = cosmo();
        let s = sigma_r(8.0, 0.0, &c).unwrap();
        assert!(
            (s - c.sigma_8).abs() < 0.01,
            "σ(8,0) = {s}, expected {}",
            c.sigma_8
        );
    }

    #[test]
    fn test_sigma_decreases_with_z() {
        let c = cosmo();
        let s0 = sigma_r(8.0, 0.0, &c).unwrap();
        let s1 = sigma_r(8.0, 1.0, &c).unwrap();
        assert!(s0 > s1, "fluctuations should decrease at higher z");
    }

    #[test]
    fn test_invalid_rejected() {
        let c = cosmo();
        assert!(transfer_function_eh(-1.0, c.omega_m, c.omega_b, c.h).is_err());
        assert!(growth_factor(-2.0, &c).is_err());
        assert!(sigma_r(0.0, 0.0, &c).is_err());
    }

    #[test]
    fn test_bao_transfer_small_k() {
        let c = cosmo();
        let t = transfer_function_eh_wiggle(0.001, c.omega_m, c.omega_b, c.h).unwrap();
        assert!(t > 0.8, "T_BAO(k→0) should approach 1, got {t}");
    }

    #[test]
    fn test_bao_transfer_large_k() {
        let c = cosmo();
        let t = transfer_function_eh_wiggle(1.0, c.omega_m, c.omega_b, c.h).unwrap();
        assert!(
            t < 0.5 && t > 0.0,
            "T_BAO(k=1) should be suppressed, got {t}"
        );
    }

    #[test]
    fn test_bao_wiggle_feature() {
        let c = cosmo();
        let mut diffs = Vec::new();
        for i in 1..50 {
            let k = 0.01 * i as f64;
            let t_nw = transfer_function_eh(k, c.omega_m, c.omega_b, c.h).unwrap();
            let t_w = transfer_function_eh_wiggle(k, c.omega_m, c.omega_b, c.h).unwrap();
            diffs.push((t_w - t_nw).abs());
        }
        let max_diff = diffs.iter().cloned().fold(0.0_f64, f64::max);
        assert!(
            max_diff > 0.001,
            "BAO wiggles should produce visible differences, max_diff={max_diff}"
        );
    }

    #[test]
    fn test_bao_transfer_bounded() {
        let c = cosmo();
        for i in 1..100 {
            let k = 0.005 * i as f64;
            let t = transfer_function_eh_wiggle(k, c.omega_m, c.omega_b, c.h).unwrap();
            assert!((0.0..=1.5).contains(&t), "T_BAO(k={k}) = {t} out of bounds");
        }
    }

    // -- dark energy --

    #[test]
    fn test_eos_lcdm() {
        let c = cosmo();
        for z in [0.0, 0.5, 1.0, 5.0] {
            let w = dark_energy_eos(z, &c).unwrap();
            assert!((w - (-1.0)).abs() < 1e-10, "ΛCDM w at z={z}: {w}");
        }
    }

    #[test]
    fn test_eos_cpl_variation() {
        let c = Cosmology::new(
            0.674, 0.315, 0.049, 0.0, 9.14e-5, 0.965, 0.811, -0.9, 0.1, 2.725,
        )
        .unwrap();
        let w0 = dark_energy_eos(0.0, &c).unwrap();
        assert!((w0 - (-0.9)).abs() < 1e-10);
        let w1 = dark_energy_eos(1.0, &c).unwrap();
        assert!(w1 > -0.9, "w should increase at z>0 for positive wa");
    }

    #[test]
    fn test_de_density_lcdm_constant() {
        let c = cosmo();
        for z in [0.5, 1.0, 3.0] {
            let r = dark_energy_density_ratio(z, &c).unwrap();
            assert!((r - 1.0).abs() < 1e-10, "ΛCDM density ratio at z={z}: {r}");
        }
    }

    #[test]
    fn test_hubble_ratio_z0() {
        let c = cosmo();
        let e = hubble_parameter_ratio(0.0, &c).unwrap();
        assert!((e - 1.0).abs() < 1e-3, "E(0) = {e}");
    }

    #[test]
    fn test_hubble_ratio_increases_with_z() {
        let c = cosmo();
        let e0 = hubble_parameter_ratio(0.0, &c).unwrap();
        let e1 = hubble_parameter_ratio(1.0, &c).unwrap();
        let e3 = hubble_parameter_ratio(3.0, &c).unwrap();
        assert!(e3 > e1 && e1 > e0);
    }

    #[test]
    fn test_growth_w_lcdm_matches_standard() {
        let c = cosmo();
        for z in [0.0, 0.5, 1.0, 2.0] {
            let d_std = growth_factor(z, &c).unwrap();
            let d_w = growth_factor_w(z, &c).unwrap();
            assert!(
                (d_std - d_w).abs() < 0.02,
                "ΛCDM growth mismatch at z={z}: std={d_std}, w={d_w}"
            );
        }
    }

    #[test]
    fn test_growth_w_phantom_faster() {
        let c_lcdm = cosmo();
        let c_phantom = Cosmology::new(
            0.674, 0.315, 0.049, 0.0, 9.14e-5, 0.965, 0.811, -1.2, 0.0, 2.725,
        )
        .unwrap();
        let d_lcdm = growth_factor_w(1.0, &c_lcdm).unwrap();
        let d_phantom = growth_factor_w(1.0, &c_phantom).unwrap();
        assert!(d_phantom > d_lcdm, "phantom DE should slow growth less");
    }

    #[test]
    fn test_comoving_distance_z0() {
        let c = cosmo();
        let chi = comoving_distance(0.0, &c).unwrap();
        assert!(chi.abs() < 1e-10);
    }

    #[test]
    fn test_comoving_distance_z1() {
        let c = cosmo();
        let chi = comoving_distance(1.0, &c).unwrap();
        assert!(chi > 2500.0 && chi < 4000.0, "χ(z=1) = {chi} Mpc");
    }

    #[test]
    fn test_comoving_distance_monotonic() {
        let c = cosmo();
        let c1 = comoving_distance(0.5, &c).unwrap();
        let c2 = comoving_distance(1.0, &c).unwrap();
        let c3 = comoving_distance(2.0, &c).unwrap();
        assert!(c3 > c2 && c2 > c1);
    }

    #[test]
    fn test_angular_power_spectrum_positive() {
        let c = cosmo();
        let cl = angular_power_spectrum_limber(10, 1.0, &c).unwrap();
        assert!(cl > 0.0, "C_l should be positive: {cl}");
    }

    #[test]
    fn test_angular_power_spectrum_decreases_with_l() {
        let c = cosmo();
        let cl10 = angular_power_spectrum_limber(10, 1.0, &c).unwrap();
        let cl1000 = angular_power_spectrum_limber(1000, 1.0, &c).unwrap();
        assert!(cl10 > cl1000, "C_10={cl10} should be > C_1000={cl1000}");
    }

    // -- distance modulus --

    #[test]
    fn test_luminosity_distance_z0() {
        let c = cosmo();
        let dl = luminosity_distance(0.0, &c).unwrap();
        assert!(dl.abs() < 1e-10);
    }

    #[test]
    fn test_luminosity_distance_increases() {
        let c = cosmo();
        let d1 = luminosity_distance(0.5, &c).unwrap();
        let d2 = luminosity_distance(1.0, &c).unwrap();
        let d3 = luminosity_distance(2.0, &c).unwrap();
        assert!(d3 > d2 && d2 > d1);
    }

    #[test]
    fn test_distance_modulus_monotonic() {
        let c = cosmo();
        let m1 = distance_modulus(0.01, &c).unwrap();
        let m2 = distance_modulus(0.1, &c).unwrap();
        let m3 = distance_modulus(1.0, &c).unwrap();
        assert!(m3 > m2 && m2 > m1, "μ should increase with z");
    }

    #[test]
    fn test_distance_modulus_z0_invalid() {
        let c = cosmo();
        assert!(distance_modulus(0.0, &c).is_err());
        assert!(distance_modulus(-1.0, &c).is_err());
    }

    // -- halofit --

    #[test]
    fn test_halofit_positive() {
        let c = cosmo();
        for k in [0.01, 0.1, 1.0, 10.0] {
            let d = halofit_nl(k, 0.0, &c).unwrap();
            assert!(d > 0.0, "Δ²_NL(k={k}) should be positive: {d}");
        }
    }

    #[test]
    fn test_halofit_increases_with_k() {
        let c = cosmo();
        let d1 = halofit_nl(0.1, 0.0, &c).unwrap();
        let d2 = halofit_nl(1.0, 0.0, &c).unwrap();
        assert!(d2 > d1, "Δ²_NL should increase with k");
    }

    #[test]
    fn test_halofit_invalid() {
        let c = cosmo();
        assert!(halofit_nl(0.0, 0.0, &c).is_err());
        assert!(halofit_nl(-1.0, 0.0, &c).is_err());
        assert!(halofit_nl(f64::NAN, 0.0, &c).is_err());
    }

    // -- sachs-wolfe --

    #[test]
    fn test_sachs_wolfe_scaling() {
        let dt = sachs_wolfe(3e-5).unwrap();
        assert!((dt - 1e-5).abs() < 1e-15);
    }

    #[test]
    fn test_isw_lcdm_nonzero() {
        let c = cosmo();
        let isw = integrated_sachs_wolfe(0.0, 2.0, &c).unwrap();
        assert!(isw.abs() > 1e-6, "ΛCDM ISW should be nonzero: {isw}");
    }

    #[test]
    fn test_isw_eds_negligible() {
        let c = Cosmology::new(
            0.674, 0.999, 0.049, 0.0, 0.0, 0.965, 0.811, -1.0, 0.0, 2.725,
        )
        .unwrap();
        let isw = integrated_sachs_wolfe(0.0, 2.0, &c).unwrap();
        assert!(isw.abs() < 0.01, "EdS ISW should be ~0: {isw}");
    }

    #[test]
    fn test_isw_invalid() {
        let c = cosmo();
        assert!(integrated_sachs_wolfe(-1.0, 2.0, &c).is_err());
        assert!(integrated_sachs_wolfe(2.0, 1.0, &c).is_err());
    }

    #[test]
    fn test_dark_energy_invalid() {
        let c = cosmo();
        assert!(dark_energy_eos(-2.0, &c).is_err());
        assert!(hubble_parameter_ratio(-2.0, &c).is_err());
        assert!(comoving_distance(-1.0, &c).is_err());
        assert!(angular_power_spectrum_limber(0, 1.0, &c).is_err());
    }

    // -- angular diameter distance --

    #[test]
    fn test_angular_diameter_distance_z0() {
        let c = cosmo();
        let da = angular_diameter_distance(0.0, &c).unwrap();
        assert!(da.abs() < 1e-10);
    }

    #[test]
    fn test_angular_diameter_distance_relation() {
        let c = cosmo();
        let z = 1.0;
        let da = angular_diameter_distance(z, &c).unwrap();
        let dl = luminosity_distance(z, &c).unwrap();
        let ratio = dl / da;
        assert!(
            (ratio - (1.0 + z).powi(2)).abs() < 1e-6,
            "d_L/d_A = {ratio}, expected {}",
            (1.0 + z).powi(2)
        );
    }

    // -- comoving volume element --

    #[test]
    fn test_comoving_volume_element_z0() {
        let c = cosmo();
        let dv = comoving_volume_element(0.0, &c).unwrap();
        assert!(dv.abs() < 1e-6, "dV/dz at z=0 should be ~0: {dv}");
    }

    #[test]
    fn test_comoving_volume_element_positive() {
        let c = cosmo();
        for z in [0.5, 1.0, 2.0, 5.0] {
            let dv = comoving_volume_element(z, &c).unwrap();
            assert!(dv > 0.0, "dV/dz should be positive at z={z}");
        }
    }

    // -- lookback time --

    #[test]
    fn test_lookback_time_z0() {
        let c = cosmo();
        let t = lookback_time(0.0, &c).unwrap();
        assert!(t.abs() < 1e-10);
    }

    #[test]
    fn test_lookback_time_monotonic() {
        let c = cosmo();
        let t1 = lookback_time(0.5, &c).unwrap();
        let t2 = lookback_time(1.0, &c).unwrap();
        let t3 = lookback_time(2.0, &c).unwrap();
        assert!(t3 > t2 && t2 > t1, "lookback time should increase with z");
    }

    #[test]
    fn test_lookback_time_z1_value() {
        let c = cosmo();
        let t = lookback_time(1.0, &c).unwrap();
        assert!(t > 6.0 && t < 10.0, "t_lb(z=1) = {t} Gyr");
    }

    // -- age of universe --

    #[test]
    fn test_age_of_universe_planck() {
        let c = cosmo();
        let age = age_of_universe(0.0, &c).unwrap();
        assert!(age > 12.0 && age < 15.0, "age = {age} Gyr, expected ~13.8");
    }

    #[test]
    fn test_age_decreases_with_z() {
        let c = cosmo();
        let a0 = age_of_universe(0.0, &c).unwrap();
        let a1 = age_of_universe(1.0, &c).unwrap();
        let a5 = age_of_universe(5.0, &c).unwrap();
        assert!(a0 > a1 && a1 > a5);
    }

    #[test]
    fn test_age_plus_lookback_equals_total() {
        let c = cosmo();
        let total = age_of_universe(0.0, &c).unwrap();
        for z in [0.5, 1.0, 2.0] {
            let age_z = age_of_universe(z, &c).unwrap();
            let tlb = lookback_time(z, &c).unwrap();
            assert!(
                (age_z + tlb - total).abs() / total < 0.02,
                "age({z}) + t_lb({z}) = {} ≠ {total}",
                age_z + tlb
            );
        }
    }

    // -- growth rate --

    #[test]
    fn test_growth_rate_z0() {
        let c = cosmo();
        let f = growth_rate(0.0, &c).unwrap();
        let expected = c.omega_m.powf(0.55);
        assert!(
            (f - expected).abs() / expected < 0.1,
            "f(z=0) = {f}, expected ~{expected}"
        );
    }

    #[test]
    fn test_growth_rate_high_z() {
        let c = cosmo();
        let f = growth_rate(10.0, &c).unwrap();
        assert!(f > 0.95, "f(z=10) = {f}, expected ~1");
    }

    #[test]
    fn test_growth_rate_increases_with_z() {
        let c = cosmo();
        let f0 = growth_rate(0.0, &c).unwrap();
        let f1 = growth_rate(1.0, &c).unwrap();
        let f5 = growth_rate(5.0, &c).unwrap();
        assert!(
            f5 > f1 && f1 > f0,
            "f should increase toward 1 at high z: f(0)={f0}, f(1)={f1}, f(5)={f5}"
        );
    }

    #[test]
    fn test_growth_rate_invalid() {
        let c = cosmo();
        assert!(growth_rate(-2.0, &c).is_err());
    }

    // -- linear power spectrum --

    #[test]
    fn test_linear_power_spectrum_positive() {
        let c = cosmo();
        let p = linear_power_spectrum(0.1, 0.0, &c).unwrap();
        assert!(p > 0.0);
    }

    #[test]
    fn test_linear_power_spectrum_decreases_with_z() {
        let c = cosmo();
        let p0 = linear_power_spectrum(0.1, 0.0, &c).unwrap();
        let p1 = linear_power_spectrum(0.1, 1.0, &c).unwrap();
        assert!(p0 > p1);
    }
}
