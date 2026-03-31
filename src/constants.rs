//! Physical constants for large-scale structure cosmology.
//!
//! All values in SI units unless noted. Sources: CODATA 2018, Planck 2018.

/// Speed of light in vacuum (m/s). CODATA 2018 exact.
pub const C: f64 = 299_792_458.0;

/// Gravitational constant G (m³ kg⁻¹ s⁻²). CODATA 2018.
pub const G: f64 = 6.674_30e-11;

/// Planck's reduced constant ℏ (J·s). CODATA 2018 exact.
pub const HBAR: f64 = 1.054_571_817e-34;

/// Boltzmann constant (J/K). CODATA 2018 exact.
pub const K_B: f64 = 1.380_649e-23;

/// Hubble constant H₀ (km/s/Mpc). Planck 2018 TT,TE,EE+lowE+lensing.
pub const H0_KM_S_MPC: f64 = 67.4;

/// 1 Megaparsec in meters. IAU 2012.
pub const MPC_M: f64 = 3.085_677_581e22;

/// 1 Megaparsec in km.
pub const MPC_KM: f64 = 3.085_677_581e19;

/// 1 kiloparsec in meters.
pub const KPC_M: f64 = 3.085_677_581e19;

/// H₀ in s⁻¹ (derived from H₀ in km/s/Mpc).
pub const H0_SI: f64 = H0_KM_S_MPC / MPC_KM;

/// Critical density today: ρ_c = 3H₀²/(8πG) (kg/m³). Derived.
pub const RHO_CRIT: f64 = 3.0 * H0_SI * H0_SI / (8.0 * std::f64::consts::PI * G);

/// Solar mass (kg). IAU 2015 nominal.
pub const M_SUN: f64 = 1.989e30;

/// Matter density parameter Ω_m. Planck 2018 TT,TE,EE+lowE+lensing.
pub const OMEGA_M: f64 = 0.315;

/// Dark energy density parameter Ω_Λ. Planck 2018 (flat ΛCDM: 1 − Ω_m).
pub const OMEGA_LAMBDA: f64 = 0.685;

/// Baryon density parameter Ω_b. Planck 2018 TT,TE,EE+lowE+lensing.
pub const OMEGA_B: f64 = 0.049;

/// Spectral index n_s. Planck 2018 TT,TE,EE+lowE+lensing.
pub const N_S: f64 = 0.965;

/// σ₈ — RMS matter fluctuations in 8 Mpc/h spheres. Planck 2018.
pub const SIGMA_8: f64 = 0.811;
