//! # Brahmanda — Large-Scale Structure Physics
//!
//! **ब्रह्माण्ड** (Sanskrit: the cosmic egg, the universe)
//!
//! Galactic structure, cosmic web topology, dark matter halos, and
//! cosmological structure formation. Part of the AGNOS science stack.
//!
//! # Modules
//!
//! - [`morphology`] — Galaxy classification (Hubble types), Sersic profiles,
//!   Faber-Jackson and Tully-Fisher scaling relations
//! - [`halo`] — Dark matter halos: NFW density profiles, virial relations,
//!   concentration-mass relations, enclosed mass and circular velocity
//! - [`cosmic_web`] — Cosmic web topology: filament/void/sheet/node classification
//!   via tidal tensor eigenvalues, density contrast, two-point correlation
//! - [`power_spectrum`] — Matter power spectrum: Eisenstein-Hu transfer function,
//!   primordial power, linear growth factor, σ(R) mass fluctuations
//!
//! # Relationship to AGNOS Science Stack
//!
//! ```text
//! hisab (math foundation)
//!   ├── hisab-mimamsa — theoretical physics (GR, QFT, cosmology, unified)
//!   │     └── unified::scale_bridge — Scale 5 consumes brahmanda outputs
//!   ├── brahmanda (this) — galactic / large-scale structure
//!   ├── tara — stellar astrophysics
//!   ├── jyotish — astronomical computation
//!   └── falak — orbital mechanics
//! ```
//!
//! # Example
//!
//! ```
//! use brahmanda::Cosmology;
//! use brahmanda::halo::HaloProperties;
//!
//! let cosmo = Cosmology::planck2018();
//! let halo = HaloProperties::from_mass(1e12, &cosmo).unwrap();
//! assert!(halo.r_vir_kpc > 150.0 && halo.r_vir_kpc < 350.0);
//! assert!(halo.concentration > 5.0 && halo.concentration < 15.0);
//! ```

pub mod constants;
pub mod cosmic_web;
pub mod cosmology;
pub mod error;
pub mod halo;
pub mod morphology;
pub mod power_spectrum;

#[cfg(feature = "logging")]
pub mod logging;

pub use cosmology::{Cosmology, FilterFunction};
pub use error::BrahmandaError;
