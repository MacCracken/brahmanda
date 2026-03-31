//! Physical constants for large-scale structure cosmology.
//!
//! All values in SI units unless noted. Sources: CODATA 2018, IAU.
//!
//! Cosmological parameters (H_0, Omega_m, etc.) are no longer here;
//! use [`crate::cosmology::Cosmology`] instead.

/// Speed of light in vacuum (m/s). CODATA 2018 exact.
pub const C: f64 = 299_792_458.0;

/// Gravitational constant G (m³ kg⁻¹ s⁻²). CODATA 2018.
pub const G: f64 = 6.674_30e-11;

/// Planck's reduced constant ℏ (J·s). CODATA 2018 exact.
pub const HBAR: f64 = 1.054_571_817e-34;

/// Boltzmann constant (J/K). CODATA 2018 exact.
pub const K_B: f64 = 1.380_649e-23;

/// 1 Megaparsec in meters. IAU 2012.
pub const MPC_M: f64 = 3.085_677_581e22;

/// 1 Megaparsec in km.
pub const MPC_KM: f64 = 3.085_677_581e19;

/// 1 kiloparsec in meters.
pub const KPC_M: f64 = 3.085_677_581e19;

/// Solar mass (kg). IAU 2015 nominal.
pub const M_SUN: f64 = 1.989e30;
