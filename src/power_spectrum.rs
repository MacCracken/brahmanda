//! Matter power spectrum — structure formation, transfer functions.
//!
//! The power spectrum P(k) describes the amplitude of density fluctuations
//! as a function of wavenumber k. Combined with the growth factor, it
//! predicts structure at all scales and redshifts.

use crate::constants::{C, H0_KM_S_MPC, N_S, OMEGA_M, SIGMA_8};
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
///
/// ```
/// use brahmanda::power_spectrum::transfer_function_eh;
///
/// // Small k (large scales): T(k) → 1
/// let t = transfer_function_eh(0.001, 0.315, 0.049, 0.674).unwrap();
/// assert!(t > 0.8);
///
/// // Large k (small scales): suppressed
/// let t = transfer_function_eh(1.0, 0.315, 0.049, 0.674).unwrap();
/// assert!(t < 0.5);
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

/// Eisenstein-Hu transfer function T(k) — full version with BAO wiggles.
///
/// Implements the baryon + CDM transfer function from Eisenstein & Hu (1998)
/// including baryon acoustic oscillations. More accurate than the no-wiggle
/// version for BAO-scale analyses.
///
/// # Arguments
/// * `k_mpc` — Wavenumber in h/Mpc.
/// * `omega_m` — Matter density parameter.
/// * `omega_b` — Baryon density parameter.
/// * `h` — Dimensionless Hubble parameter (H₀/100).
///
/// ```
/// use brahmanda::power_spectrum::transfer_function_eh_wiggle;
///
/// let t = transfer_function_eh_wiggle(0.001, 0.315, 0.049, 0.674).unwrap();
/// assert!(t > 0.8); // large scales → T ≈ 1
///
/// // Should show oscillatory features relative to no-wiggle
/// let t_bao = transfer_function_eh_wiggle(0.1, 0.315, 0.049, 0.674).unwrap();
/// assert!(t_bao > 0.0 && t_bao < 1.0);
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
    let k_silk = 1.6 * omega_bh2.powf(0.52) * omega_mh2.powf(0.73)
        * (1.0 + (10.4 * omega_mh2).powf(-0.95));

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
    let g_y = y
        * (-6.0 * sqrt_1py
            + (2.0 + 3.0 * y) * ((sqrt_1py + 1.0) / (sqrt_1py - 1.0)).ln());

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
        + alpha_b / (1.0 + (beta_b / ks).powi(3))
            * (-(k_mpc / k_silk).powf(1.4)).exp())
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
/// use brahmanda::power_spectrum::primordial_power;
///
/// let p = primordial_power(1.0, 0.965).unwrap();
/// assert!((p - 1.0).abs() < 1e-10); // k=1 → k^n_s = 1
/// ```
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
///
/// ```
/// use brahmanda::power_spectrum::growth_factor;
///
/// // D(z=0) = 1 by definition
/// let d = growth_factor(0.0, 0.315).unwrap();
/// assert!((d - 1.0).abs() < 1e-10);
///
/// // Structure grows: D decreases at higher z
/// let d1 = growth_factor(1.0, 0.315).unwrap();
/// assert!(d1 < d);
/// ```
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
///
/// ```
/// use brahmanda::power_spectrum::sigma_r;
///
/// // σ(R=8, z=0) ≈ σ₈
/// let s = sigma_r(8.0, 0.0).unwrap();
/// assert!((s - 0.811).abs() < 0.01);
///
/// // Fluctuations decrease at higher z
/// let s1 = sigma_r(8.0, 1.0).unwrap();
/// assert!(s1 < s);
/// ```
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

/// Dark energy equation of state — CPL parameterization.
///
/// w(z) = w₀ + w_a × z/(1+z)
///
/// Chevallier-Polarski-Linder (CPL) parameterization. ΛCDM corresponds
/// to w₀ = -1, w_a = 0.
///
/// # Arguments
/// * `z` — Redshift.
/// * `w0` — Equation of state at z=0 (w₀ = -1 for ΛCDM).
/// * `wa` — Time evolution parameter (w_a = 0 for ΛCDM).
///
/// ```
/// use brahmanda::power_spectrum::dark_energy_eos;
///
/// // ΛCDM: w = -1 at all redshifts
/// let w = dark_energy_eos(0.5, -1.0, 0.0).unwrap();
/// assert!((w - (-1.0)).abs() < 1e-10);
///
/// // Thawing model: w > -1 at high z
/// let w = dark_energy_eos(1.0, -0.9, 0.1).unwrap();
/// assert!(w > -1.0);
/// ```
#[inline]
pub fn dark_energy_eos(z: f64, w0: f64, wa: f64) -> Result<f64, BrahmandaError> {
    require_finite(z, "dark_energy_eos")?;
    require_finite(w0, "dark_energy_eos")?;
    require_finite(wa, "dark_energy_eos")?;
    if z < -1.0 {
        return Err(BrahmandaError::InvalidStructure(
            "dark_energy_eos: z must be >= -1".to_string(),
        ));
    }
    ensure_finite(w0 + wa * z / (1.0 + z), "dark_energy_eos")
}

/// Dark energy density evolution factor for CPL parameterization.
///
/// ρ_DE(z)/ρ_DE(0) = (1+z)^(3(1+w₀+w_a)) × exp(-3 w_a z/(1+z))
///
/// ```
/// use brahmanda::power_spectrum::dark_energy_density_ratio;
///
/// // ΛCDM: constant dark energy density
/// let ratio = dark_energy_density_ratio(1.0, -1.0, 0.0).unwrap();
/// assert!((ratio - 1.0).abs() < 1e-10);
/// ```
pub fn dark_energy_density_ratio(z: f64, w0: f64, wa: f64) -> Result<f64, BrahmandaError> {
    require_finite(z, "dark_energy_density_ratio")?;
    require_finite(w0, "dark_energy_density_ratio")?;
    require_finite(wa, "dark_energy_density_ratio")?;
    if z < -1.0 {
        return Err(BrahmandaError::InvalidStructure(
            "dark_energy_density_ratio: z must be >= -1".to_string(),
        ));
    }
    let a = 1.0 / (1.0 + z);
    let exponent = 3.0 * (1.0 + w0 + wa);
    let ratio = (1.0 + z).powf(exponent) * (-3.0 * wa * (1.0 - a)).exp();
    ensure_finite(ratio, "dark_energy_density_ratio")
}

/// Hubble parameter E(z) = H(z)/H₀ for w₀w_a dark energy.
///
/// E²(z) = Ω_m(1+z)³ + Ω_DE × (ρ_DE(z)/ρ_DE(0))
///
/// ```
/// use brahmanda::power_spectrum::hubble_parameter_ratio;
///
/// // At z=0, E(0) = 1 by definition
/// let e = hubble_parameter_ratio(0.0, 0.315, -1.0, 0.0).unwrap();
/// assert!((e - 1.0).abs() < 1e-6);
/// ```
pub fn hubble_parameter_ratio(
    z: f64,
    omega_m: f64,
    w0: f64,
    wa: f64,
) -> Result<f64, BrahmandaError> {
    require_finite(z, "hubble_parameter_ratio")?;
    require_finite(omega_m, "hubble_parameter_ratio")?;
    if z < -1.0 {
        return Err(BrahmandaError::InvalidStructure(
            "hubble_parameter_ratio: z must be >= -1".to_string(),
        ));
    }
    let omega_de = 1.0 - omega_m;
    let de_ratio = dark_energy_density_ratio(z, w0, wa)?;
    let e2 = omega_m * (1.0 + z).powi(3) + omega_de * de_ratio;
    ensure_finite(e2.sqrt(), "hubble_parameter_ratio")
}

/// Linear growth factor for w₀w_a dark energy — CPT92 fitting formula.
///
/// Generalizes the Carroll, Press & Turner (1992) approximation to
/// w₀w_a dark energy by computing Ω_m(z) and Ω_DE(z) at each redshift.
///
/// Returns D(z)/D(0).
///
/// ```
/// use brahmanda::power_spectrum::growth_factor_w;
///
/// // ΛCDM should match the standard growth factor
/// let d = growth_factor_w(1.0, 0.315, -1.0, 0.0).unwrap();
/// assert!(d > 0.0 && d < 1.0);
///
/// // Phantom dark energy (w < -1) slows growth less
/// let d_phantom = growth_factor_w(1.0, 0.315, -1.2, 0.0).unwrap();
/// assert!(d_phantom > d);
/// ```
pub fn growth_factor_w(
    z: f64,
    omega_m: f64,
    w0: f64,
    wa: f64,
) -> Result<f64, BrahmandaError> {
    require_finite(z, "growth_factor_w")?;
    require_finite(omega_m, "growth_factor_w")?;
    require_finite(w0, "growth_factor_w")?;
    require_finite(wa, "growth_factor_w")?;
    if z < -1.0 {
        return Err(BrahmandaError::InvalidStructure(
            "growth_factor_w: z must be >= -1".to_string(),
        ));
    }

    // Use the CPT-style approximation generalized for w(z):
    // At each z, compute Ω_m(z) and Ω_DE(z), then use the fitting formula
    let growth_at = |zz: f64| -> Result<f64, BrahmandaError> {
        let a = 1.0 / (1.0 + zz);
        let e = hubble_parameter_ratio(zz, omega_m, w0, wa)?;
        let e2 = e * e;
        let om_z = omega_m * (1.0 + zz).powi(3) / e2;
        let ol_z = 1.0 - om_z;
        // CPT92 approximation (works well for slowly-varying w)
        let d = 2.5 * om_z
            / (om_z.powf(4.0 / 7.0) - ol_z + (1.0 + om_z / 2.0) * (1.0 + ol_z / 70.0));
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
/// Uses Simpson's rule numerical integration.
///
/// ```
/// use brahmanda::power_spectrum::comoving_distance;
///
/// let chi = comoving_distance(0.0, 0.315, -1.0, 0.0).unwrap();
/// assert!(chi.abs() < 1e-10); // zero distance at z=0
///
/// let chi = comoving_distance(1.0, 0.315, -1.0, 0.0).unwrap();
/// assert!(chi > 2000.0 && chi < 5000.0); // ~3300 Mpc
/// ```
pub fn comoving_distance(
    z: f64,
    omega_m: f64,
    w0: f64,
    wa: f64,
) -> Result<f64, BrahmandaError> {
    require_finite(z, "comoving_distance")?;
    if z < 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "comoving_distance: z must be non-negative".to_string(),
        ));
    }
    if z == 0.0 {
        return Ok(0.0);
    }

    // Simpson's rule with N steps
    let n = 200_usize;
    let dz = z / n as f64;
    let mut sum = 0.0;

    for i in 0..n {
        let z0 = i as f64 * dz;
        let z1 = z0 + dz / 2.0;
        let z2 = z0 + dz;
        let e0 = hubble_parameter_ratio(z0, omega_m, w0, wa)?;
        let e1 = hubble_parameter_ratio(z1, omega_m, w0, wa)?;
        let e2 = hubble_parameter_ratio(z2, omega_m, w0, wa)?;
        sum += (dz / 6.0) * (1.0 / e0 + 4.0 / e1 + 1.0 / e2);
    }

    // c/H₀ in Mpc: c in km/s, H₀ in km/s/Mpc → c/H₀ in Mpc
    let c_over_h0 = (C / 1e3) / H0_KM_S_MPC;
    ensure_finite(c_over_h0 * sum, "comoving_distance")
}

/// Angular power spectrum C_l — Limber approximation.
///
/// C_l = ∫ dχ [W(χ)]² P(k=l/χ, z(χ)) / χ²
///
/// Uses a simplified computation with the matter power spectrum P(k) ∝ k^n_s T²(k)
/// and a galaxy/matter window function. The result is the unnormalized
/// angular power spectrum for matter density projections.
///
/// # Arguments
/// * `l` — Multipole moment (must be > 0).
/// * `z_max` — Maximum redshift of the source distribution.
/// * `omega_m` — Matter density parameter.
/// * `omega_b` — Baryon density parameter.
/// * `h` — Dimensionless Hubble parameter.
/// * `w0` — Dark energy equation of state w₀.
/// * `wa` — Dark energy evolution parameter w_a.
///
/// ```
/// use brahmanda::power_spectrum::angular_power_spectrum_limber;
///
/// let cl_10 = angular_power_spectrum_limber(10, 1.0, 0.315, 0.049, 0.674, -1.0, 0.0).unwrap();
/// let cl_100 = angular_power_spectrum_limber(100, 1.0, 0.315, 0.049, 0.674, -1.0, 0.0).unwrap();
/// assert!(cl_10 > 0.0);
/// assert!(cl_100 > 0.0);
/// ```
pub fn angular_power_spectrum_limber(
    l: u32,
    z_max: f64,
    omega_m: f64,
    omega_b: f64,
    h: f64,
    w0: f64,
    wa: f64,
) -> Result<f64, BrahmandaError> {
    require_finite(z_max, "angular_power_spectrum_limber")?;
    require_finite(omega_m, "angular_power_spectrum_limber")?;
    require_finite(omega_b, "angular_power_spectrum_limber")?;
    require_finite(h, "angular_power_spectrum_limber")?;
    if l == 0 {
        return Err(BrahmandaError::InvalidStructure(
            "angular_power_spectrum_limber: l must be > 0".to_string(),
        ));
    }
    if z_max <= 0.0 || omega_m <= 0.0 || h <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "angular_power_spectrum_limber: z_max, omega_m, h must be positive".to_string(),
        ));
    }

    let l_f = l as f64;
    let chi_max = comoving_distance(z_max, omega_m, w0, wa)?;

    // Uniform window function W(χ) = 1/χ_max for a flat source distribution
    // Then C_l = (1/χ_max²) ∫₀^χ_max dχ P(k=l/χ) / χ²

    // Integrate in redshift space using Simpson's rule
    let n = 100_usize;
    let dz = z_max / n as f64;
    let mut sum = 0.0;

    let integrand = |z: f64| -> Result<f64, BrahmandaError> {
        if z < 1e-6 {
            return Ok(0.0); // skip z=0 singularity
        }
        let chi = comoving_distance(z, omega_m, w0, wa)?;
        if chi < 1e-6 {
            return Ok(0.0);
        }
        let k = (l_f + 0.5) / chi; // Limber: k = (l+1/2)/χ in Mpc⁻¹
        let k_h = k / h; // convert physical Mpc⁻¹ → h/Mpc

        if k_h <= 0.0 || k_h > 100.0 {
            return Ok(0.0);
        }

        let t = transfer_function_eh(k_h, omega_m, omega_b, h)?;
        let d = growth_factor_w(z, omega_m, w0, wa)?;

        // P(k,z) ∝ k^n_s × T²(k) × D²(z)
        let pk = k_h.powf(N_S) * t * t * d * d;

        let e = hubble_parameter_ratio(z, omega_m, w0, wa)?;

        // dχ/dz = c/(H₀ E(z)), and integrand is P(k)/χ² × (dχ/dz)
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

    let c_over_h0 = (C / 1e3) / H0_KM_S_MPC;
    let cl = c_over_h0 * sum / (chi_max * chi_max);
    ensure_finite(cl, "angular_power_spectrum_limber")
}

/// Luminosity distance d_L(z) in Mpc.
///
/// d_L = (1+z) × χ(z) for a flat universe.
///
/// ```
/// use brahmanda::power_spectrum::luminosity_distance;
///
/// let dl = luminosity_distance(0.0, 0.315, -1.0, 0.0).unwrap();
/// assert!(dl.abs() < 1e-10);
///
/// let dl = luminosity_distance(1.0, 0.315, -1.0, 0.0).unwrap();
/// assert!(dl > 5000.0 && dl < 8000.0); // ~6600 Mpc
/// ```
pub fn luminosity_distance(
    z: f64,
    omega_m: f64,
    w0: f64,
    wa: f64,
) -> Result<f64, BrahmandaError> {
    let chi = comoving_distance(z, omega_m, w0, wa)?;
    ensure_finite((1.0 + z) * chi, "luminosity_distance")
}

/// Distance modulus μ(z) for SN Ia cosmology.
///
/// μ = 5 log₁₀(d_L / 10 pc) = 5 log₁₀(d_L [Mpc]) + 25
///
/// ```
/// use brahmanda::power_spectrum::distance_modulus;
///
/// // At z≈0.01, μ ≈ 33 (nearby SN Ia)
/// let mu = distance_modulus(0.01, 0.315, -1.0, 0.0).unwrap();
/// assert!(mu > 30.0 && mu < 38.0, "μ(z=0.01) = {mu}");
///
/// // At z=1, μ ≈ 44
/// let mu = distance_modulus(1.0, 0.315, -1.0, 0.0).unwrap();
/// assert!(mu > 42.0 && mu < 46.0, "μ(z=1) = {mu}");
/// ```
pub fn distance_modulus(
    z: f64,
    omega_m: f64,
    w0: f64,
    wa: f64,
) -> Result<f64, BrahmandaError> {
    require_finite(z, "distance_modulus")?;
    if z <= 0.0 {
        return Err(BrahmandaError::InvalidStructure(
            "distance_modulus: z must be positive".to_string(),
        ));
    }
    let dl = luminosity_distance(z, omega_m, w0, wa)?;
    // μ = 5 log₁₀(d_L [pc] / 10) = 5 log₁₀(d_L [Mpc] × 10⁶) + 25 - 5
    // Simplified: μ = 5 log₁₀(d_L [Mpc]) + 25
    ensure_finite(5.0 * dl.log10() + 25.0, "distance_modulus")
}

/// Effective spectral index n_eff at scale R (Mpc/h).
///
/// n_eff = -2 × d ln σ(R) / d ln R - 3
///
/// Computed via finite differencing of σ(R).
fn effective_spectral_index(r_mpc: f64, z: f64) -> Result<f64, BrahmandaError> {
    let dr = r_mpc * 0.01;
    let s_lo = sigma_r(r_mpc - dr, z)?;
    let s_hi = sigma_r(r_mpc + dr, z)?;
    let dlns_dlnr = (s_hi.ln() - s_lo.ln()) / ((r_mpc + dr).ln() - (r_mpc - dr).ln());
    ensure_finite(-2.0 * dlns_dlnr - 3.0, "effective_spectral_index")
}

/// Nonlinear power spectrum — Smith et al. (2003) Halofit.
///
/// Maps the linear dimensionless power Δ²_lin(k) to the nonlinear Δ²_NL(k)
/// using the Halofit fitting formula.
///
/// Returns the nonlinear dimensionless power spectrum Δ²_NL(k) = k³ P_NL(k) / (2π²).
///
/// # Arguments
/// * `k_mpc` — Wavenumber in h/Mpc.
/// * `z` — Redshift.
///
/// ```
/// use brahmanda::power_spectrum::halofit_nl;
///
/// // Nonlinear power is positive and finite
/// let dnl = halofit_nl(0.1, 0.0).unwrap();
/// assert!(dnl > 0.0);
///
/// // Small scales have more nonlinear power than large scales
/// let dnl_small = halofit_nl(1.0, 0.0).unwrap();
/// assert!(dnl_small > dnl);
/// ```
pub fn halofit_nl(k_mpc: f64, z: f64) -> Result<f64, BrahmandaError> {
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
        let s = sigma_r(r_mid, z)?;
        if s > 1.0 {
            r_lo = r_mid;
        } else {
            r_hi = r_mid;
        }
    }
    let r_nl = (r_lo * r_hi).sqrt();
    let k_sigma = 1.0 / r_nl;

    let n_eff = effective_spectral_index(r_nl, z)?;
    let n2 = n_eff * n_eff;

    let y = k_mpc / k_sigma;

    // Smith et al. 2003, Table 2
    let a_n = 10.0_f64.powf(
        1.4861 + 1.8369 * n_eff + 1.6762 * n2 + 0.7940 * n_eff * n2
            + 0.1670 * n2 * n2 - 0.6206 * n_eff.powi(5),
    );
    let b_n = 10.0_f64.powf(0.9463 + 0.9466 * n_eff + 0.3084 * n2 - 0.9400 * n_eff * n2);
    let c_n = 10.0_f64.powf(-0.2807 + 0.6669 * n_eff + 0.3214 * n2 - 0.0793 * n_eff * n2);
    let mu_n = 10.0_f64.powf(-3.5442 + 0.1908 * n_eff);
    let nu_n = 10.0_f64.powf(0.9589 + 1.2857 * n_eff);
    let alpha_n = (1.3884 + 0.3700 * n_eff - 0.1452 * n2).abs();
    let beta_n = 0.8291 + 0.9854 * n_eff + 0.3401 * n2;

    // Linear Δ²(k) ∝ σ²(1/k) approximately
    let sigma_k = sigma_r(1.0 / k_mpc, z)?;
    let delta2_lin = sigma_k * sigma_k;

    // Two-halo (quasi-linear) term: Δ²_Q
    let f_y = y / 4.0 + y * y / 8.0;
    let delta2_q = delta2_lin
        * ((1.0 + delta2_lin).powf(beta_n) / (1.0 + alpha_n * delta2_lin))
        * (-f_y).exp();

    // One-halo term: Δ²_H
    let delta2_h = a_n * y.powf(3.0 * mu_n) / (1.0 + b_n * y.powf(-1.0))
        / (1.0 + c_n * y.powf(3.0 * nu_n));

    ensure_finite(delta2_q + delta2_h, "halofit_nl")
}

/// Ordinary Sachs-Wolfe effect: ΔT/T from gravitational potential.
///
/// (ΔT/T)_SW = Φ / 3
///
/// where Φ is the Newtonian gravitational potential at last scattering.
/// Valid for large angular scales (l < 30).
///
/// # Arguments
/// * `phi` — Gravitational potential (dimensionless, Φ/c²).
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
/// (ΔT/T)_ISW = -2 ∫ (dΦ/dt) dt/c = 2 ∫ (dΦ/dz) dz
///
/// In ΛCDM, the ISW arises because Φ decays at late times when dark
/// energy dominates. Uses the growth rate f = dlnD/dlna:
///
/// ISW ∝ ∫ (f - 1) D(z) H(z)/H₀ dz
///
/// In EdS (Ω_m = 1), f = 1 so ISW vanishes identically.
///
/// Returns the unnormalized ISW contribution between z_min and z_max.
///
/// ```
/// use brahmanda::power_spectrum::integrated_sachs_wolfe;
///
/// // ISW is non-zero in ΛCDM (Ω_Λ > 0)
/// let isw = integrated_sachs_wolfe(0.0, 2.0, 0.315, -1.0, 0.0).unwrap();
/// assert!(isw.abs() > 0.0);
///
/// // EdS universe (Ω_m=1): no ISW effect (potentials constant)
/// let isw_eds = integrated_sachs_wolfe(0.0, 2.0, 1.0, -1.0, 0.0).unwrap();
/// assert!(isw_eds.abs() < 0.01, "EdS ISW should be ~0: {isw_eds}");
/// ```
pub fn integrated_sachs_wolfe(
    z_min: f64,
    z_max: f64,
    omega_m: f64,
    w0: f64,
    wa: f64,
) -> Result<f64, BrahmandaError> {
    require_finite(z_min, "integrated_sachs_wolfe")?;
    require_finite(z_max, "integrated_sachs_wolfe")?;
    if z_min < 0.0 || z_max <= z_min {
        return Err(BrahmandaError::InvalidStructure(
            "integrated_sachs_wolfe: need 0 <= z_min < z_max".to_string(),
        ));
    }

    // ISW ∝ ∫ (f(z) - 1) × D(z) × H(z)/H₀ dz
    // where f = dlnD/dlna is the growth rate.
    // In EdS, f = 1 exactly so ISW = 0. In ΛCDM, f < 1 at late times.
    let n = 200_usize;
    let dz = (z_max - z_min) / n as f64;
    let eps = 1e-4;
    let mut sum = 0.0;

    let isw_integrand = |z: f64| -> Result<f64, BrahmandaError> {
        let d = growth_factor_w(z, omega_m, w0, wa)?;
        let e = hubble_parameter_ratio(z, omega_m, w0, wa)?;

        // f = dlnD/dlna = -(1+z) dlnD/dz
        let d_hi = growth_factor_w(z + eps, omega_m, w0, wa)?;
        let d_lo = growth_factor_w((z - eps).max(0.0), omega_m, w0, wa)?;
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
    use crate::constants::OMEGA_B;

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

    #[test]
    fn test_bao_transfer_small_k() {
        let t = transfer_function_eh_wiggle(0.001, OMEGA_M, OMEGA_B, 0.674).unwrap();
        assert!(t > 0.8, "T_BAO(k→0) should approach 1, got {t}");
    }

    #[test]
    fn test_bao_transfer_large_k() {
        let t = transfer_function_eh_wiggle(1.0, OMEGA_M, OMEGA_B, 0.674).unwrap();
        assert!(t < 0.5 && t > 0.0, "T_BAO(k=1) should be suppressed, got {t}");
    }

    #[test]
    fn test_bao_wiggle_feature() {
        // The BAO version should differ from the no-wiggle at intermediate k
        // due to oscillatory features around k ~ 0.05-0.3 h/Mpc
        let mut diffs = Vec::new();
        for i in 1..50 {
            let k = 0.01 * i as f64;
            let t_nw = transfer_function_eh(k, OMEGA_M, OMEGA_B, 0.674).unwrap();
            let t_w = transfer_function_eh_wiggle(k, OMEGA_M, OMEGA_B, 0.674).unwrap();
            diffs.push((t_w - t_nw).abs());
        }
        let max_diff = diffs.iter().cloned().fold(0.0_f64, f64::max);
        assert!(max_diff > 0.001, "BAO wiggles should produce visible differences, max_diff={max_diff}");
    }

    #[test]
    fn test_bao_transfer_bounded() {
        for i in 1..100 {
            let k = 0.005 * i as f64;
            let t = transfer_function_eh_wiggle(k, OMEGA_M, OMEGA_B, 0.674).unwrap();
            assert!((0.0..=1.5).contains(&t), "T_BAO(k={k}) = {t} out of bounds");
        }
    }

    // -- dark energy --

    #[test]
    fn test_eos_lcdm() {
        for z in [0.0, 0.5, 1.0, 5.0] {
            let w = dark_energy_eos(z, -1.0, 0.0).unwrap();
            assert!((w - (-1.0)).abs() < 1e-10, "ΛCDM w at z={z}: {w}");
        }
    }

    #[test]
    fn test_eos_cpl_variation() {
        let w0 = dark_energy_eos(0.0, -0.9, 0.1).unwrap();
        assert!((w0 - (-0.9)).abs() < 1e-10);
        let w1 = dark_energy_eos(1.0, -0.9, 0.1).unwrap();
        assert!(w1 > -0.9, "w should increase at z>0 for positive wa");
    }

    #[test]
    fn test_de_density_lcdm_constant() {
        for z in [0.5, 1.0, 3.0] {
            let r = dark_energy_density_ratio(z, -1.0, 0.0).unwrap();
            assert!((r - 1.0).abs() < 1e-10, "ΛCDM density ratio at z={z}: {r}");
        }
    }

    #[test]
    fn test_hubble_ratio_z0() {
        let e = hubble_parameter_ratio(0.0, OMEGA_M, -1.0, 0.0).unwrap();
        assert!((e - 1.0).abs() < 1e-6, "E(0) = {e}");
    }

    #[test]
    fn test_hubble_ratio_increases_with_z() {
        let e0 = hubble_parameter_ratio(0.0, OMEGA_M, -1.0, 0.0).unwrap();
        let e1 = hubble_parameter_ratio(1.0, OMEGA_M, -1.0, 0.0).unwrap();
        let e3 = hubble_parameter_ratio(3.0, OMEGA_M, -1.0, 0.0).unwrap();
        assert!(e3 > e1 && e1 > e0);
    }

    #[test]
    fn test_growth_w_lcdm_matches_standard() {
        for z in [0.0, 0.5, 1.0, 2.0] {
            let d_std = growth_factor(z, OMEGA_M).unwrap();
            let d_w = growth_factor_w(z, OMEGA_M, -1.0, 0.0).unwrap();
            assert!(
                (d_std - d_w).abs() < 0.02,
                "ΛCDM growth mismatch at z={z}: std={d_std}, w={d_w}"
            );
        }
    }

    #[test]
    fn test_growth_w_phantom_faster() {
        let d_lcdm = growth_factor_w(1.0, OMEGA_M, -1.0, 0.0).unwrap();
        let d_phantom = growth_factor_w(1.0, OMEGA_M, -1.2, 0.0).unwrap();
        assert!(d_phantom > d_lcdm, "phantom DE should slow growth less");
    }

    #[test]
    fn test_comoving_distance_z0() {
        let chi = comoving_distance(0.0, OMEGA_M, -1.0, 0.0).unwrap();
        assert!(chi.abs() < 1e-10);
    }

    #[test]
    fn test_comoving_distance_z1() {
        let chi = comoving_distance(1.0, OMEGA_M, -1.0, 0.0).unwrap();
        // For Planck 2018 ΛCDM, χ(z=1) ≈ 3300 Mpc
        assert!(chi > 2500.0 && chi < 4000.0, "χ(z=1) = {chi} Mpc");
    }

    #[test]
    fn test_comoving_distance_monotonic() {
        let c1 = comoving_distance(0.5, OMEGA_M, -1.0, 0.0).unwrap();
        let c2 = comoving_distance(1.0, OMEGA_M, -1.0, 0.0).unwrap();
        let c3 = comoving_distance(2.0, OMEGA_M, -1.0, 0.0).unwrap();
        assert!(c3 > c2 && c2 > c1);
    }

    #[test]
    fn test_angular_power_spectrum_positive() {
        let cl = angular_power_spectrum_limber(10, 1.0, OMEGA_M, OMEGA_B, 0.674, -1.0, 0.0)
            .unwrap();
        assert!(cl > 0.0, "C_l should be positive: {cl}");
    }

    #[test]
    fn test_angular_power_spectrum_decreases_with_l() {
        let cl10 = angular_power_spectrum_limber(10, 1.0, OMEGA_M, OMEGA_B, 0.674, -1.0, 0.0)
            .unwrap();
        let cl1000 = angular_power_spectrum_limber(1000, 1.0, OMEGA_M, OMEGA_B, 0.674, -1.0, 0.0)
            .unwrap();
        // At high l, C_l should generally decrease (modulo BAO features)
        assert!(cl10 > cl1000, "C_10={cl10} should be > C_1000={cl1000}");
    }

    // -- distance modulus --

    #[test]
    fn test_luminosity_distance_z0() {
        let dl = luminosity_distance(0.0, OMEGA_M, -1.0, 0.0).unwrap();
        assert!(dl.abs() < 1e-10);
    }

    #[test]
    fn test_luminosity_distance_increases() {
        let d1 = luminosity_distance(0.5, OMEGA_M, -1.0, 0.0).unwrap();
        let d2 = luminosity_distance(1.0, OMEGA_M, -1.0, 0.0).unwrap();
        let d3 = luminosity_distance(2.0, OMEGA_M, -1.0, 0.0).unwrap();
        assert!(d3 > d2 && d2 > d1);
    }

    #[test]
    fn test_distance_modulus_monotonic() {
        let m1 = distance_modulus(0.01, OMEGA_M, -1.0, 0.0).unwrap();
        let m2 = distance_modulus(0.1, OMEGA_M, -1.0, 0.0).unwrap();
        let m3 = distance_modulus(1.0, OMEGA_M, -1.0, 0.0).unwrap();
        assert!(m3 > m2 && m2 > m1, "μ should increase with z");
    }

    #[test]
    fn test_distance_modulus_z0_invalid() {
        assert!(distance_modulus(0.0, OMEGA_M, -1.0, 0.0).is_err());
        assert!(distance_modulus(-1.0, OMEGA_M, -1.0, 0.0).is_err());
    }

    // -- halofit --

    #[test]
    fn test_halofit_positive() {
        for k in [0.01, 0.1, 1.0, 10.0] {
            let d = halofit_nl(k, 0.0).unwrap();
            assert!(d > 0.0, "Δ²_NL(k={k}) should be positive: {d}");
        }
    }

    #[test]
    fn test_halofit_increases_with_k() {
        let d1 = halofit_nl(0.1, 0.0).unwrap();
        let d2 = halofit_nl(1.0, 0.0).unwrap();
        assert!(d2 > d1, "Δ²_NL should increase with k");
    }

    #[test]
    fn test_halofit_invalid() {
        assert!(halofit_nl(0.0, 0.0).is_err());
        assert!(halofit_nl(-1.0, 0.0).is_err());
        assert!(halofit_nl(f64::NAN, 0.0).is_err());
    }

    // -- sachs-wolfe --

    #[test]
    fn test_sachs_wolfe_scaling() {
        let dt = sachs_wolfe(3e-5).unwrap();
        assert!((dt - 1e-5).abs() < 1e-15);
    }

    #[test]
    fn test_isw_lcdm_nonzero() {
        let isw = integrated_sachs_wolfe(0.0, 2.0, OMEGA_M, -1.0, 0.0).unwrap();
        assert!(isw.abs() > 1e-6, "ΛCDM ISW should be nonzero: {isw}");
    }

    #[test]
    fn test_isw_eds_negligible() {
        // In EdS (Ω_m=1), potentials are constant → ISW ≈ 0
        let isw = integrated_sachs_wolfe(0.0, 2.0, 0.999, -1.0, 0.0).unwrap();
        assert!(isw.abs() < 0.01, "EdS ISW should be ~0: {isw}");
    }

    #[test]
    fn test_isw_invalid() {
        assert!(integrated_sachs_wolfe(-1.0, 2.0, OMEGA_M, -1.0, 0.0).is_err());
        assert!(integrated_sachs_wolfe(2.0, 1.0, OMEGA_M, -1.0, 0.0).is_err());
    }

    #[test]
    fn test_dark_energy_invalid() {
        assert!(dark_energy_eos(-2.0, -1.0, 0.0).is_err());
        assert!(dark_energy_eos(f64::NAN, -1.0, 0.0).is_err());
        assert!(hubble_parameter_ratio(-2.0, OMEGA_M, -1.0, 0.0).is_err());
        assert!(comoving_distance(-1.0, OMEGA_M, -1.0, 0.0).is_err());
        assert!(angular_power_spectrum_limber(0, 1.0, OMEGA_M, OMEGA_B, 0.674, -1.0, 0.0).is_err());
    }
}
