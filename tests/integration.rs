//! Cross-module integration tests.

use brahmanda::Cosmology;
use brahmanda::cosmic_web::{self, HessianMorphology, WebEnvironment};
use brahmanda::halo;
use brahmanda::morphology::{self, HubbleType};
use brahmanda::power_spectrum;

fn cosmo() -> Cosmology {
    Cosmology::planck2018()
}

#[test]
fn test_milky_way_halo_and_galaxy() {
    let c = cosmo();
    let halo = halo::HaloProperties::from_mass(1e12, &c).unwrap();
    assert!(halo.r_vir_kpc > 150.0 && halo.r_vir_kpc < 350.0);

    let rho_s = 1e7;
    let v = halo::nfw_circular_velocity(8.0, rho_s, halo.r_s_kpc).unwrap();
    assert!(v > 50.0 && v < 500.0, "v_c at 8 kpc: {v} km/s");
}

#[test]
fn test_structure_formation_consistency() {
    let c = cosmo();
    let sigma_0 = power_spectrum::sigma_r(8.0, 0.0, &c).unwrap();
    let sigma_2 = power_spectrum::sigma_r(8.0, 2.0, &c).unwrap();
    assert!(sigma_0 > sigma_2, "less structure at high z");
}

#[test]
fn test_cosmic_web_overdense_is_node() {
    let env = cosmic_web::classify_web_environment(&[2.0, 1.5, 0.5], 0.0).unwrap();
    assert_eq!(env, cosmic_web::WebEnvironment::Node);

    let delta = cosmic_web::density_contrast(10.0, 1.0).unwrap();
    assert!(delta > 0.0, "nodes are overdense");
}

#[test]
fn test_sersic_elliptical_vs_spiral() {
    let i_elliptical = morphology::sersic_profile(0.1, 1.0, 4.0).unwrap();
    let i_spiral = morphology::sersic_profile(0.1, 1.0, 1.0).unwrap();
    assert!(i_elliptical > i_spiral, "elliptical brighter in center");
}

// -- Display trait tests --

#[test]
fn hubble_type_display() {
    assert_eq!(HubbleType::Elliptical.to_string(), "Elliptical");
    assert_eq!(HubbleType::BarredSpiral.to_string(), "Barred Spiral");
    assert_eq!(HubbleType::Irregular.to_string(), "Irregular");
}

#[test]
fn web_environment_display() {
    assert_eq!(WebEnvironment::Node.to_string(), "Node");
    assert_eq!(WebEnvironment::Void.to_string(), "Void");
}

#[test]
fn hessian_morphology_display() {
    assert_eq!(HessianMorphology::Blob.to_string(), "Blob");
    assert_eq!(HessianMorphology::Filament.to_string(), "Filament");
    assert_eq!(HessianMorphology::Wall.to_string(), "Wall");
    assert_eq!(HessianMorphology::Background.to_string(), "Background");
}
