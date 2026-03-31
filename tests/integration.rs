//! Cross-module integration tests.

use brahmanda::cosmic_web;
use brahmanda::halo;
use brahmanda::morphology;
use brahmanda::power_spectrum;

#[test]
fn test_milky_way_halo_and_galaxy() {
    // MW halo: ~10¹² M_sun
    let halo = halo::HaloProperties::from_mass(1e12).unwrap();
    assert!(halo.r_vir_kpc > 150.0 && halo.r_vir_kpc < 350.0);

    // NFW velocity curve should be reasonable at ~8 kpc (solar radius)
    let rho_s = 1e7; // typical characteristic density
    let v = halo::nfw_circular_velocity(8.0, rho_s, halo.r_s_kpc).unwrap();
    assert!(v > 50.0 && v < 500.0, "v_c at 8 kpc: {v} km/s");
}

#[test]
fn test_structure_formation_consistency() {
    // Growth factor decreases at higher z → σ decreases → fewer structures
    let sigma_0 = power_spectrum::sigma_r(8.0, 0.0).unwrap();
    let sigma_2 = power_spectrum::sigma_r(8.0, 2.0).unwrap();
    assert!(sigma_0 > sigma_2, "less structure at high z");
}

#[test]
fn test_cosmic_web_overdense_is_node() {
    // In a node: all eigenvalues positive → density contrast >> 0
    let env = cosmic_web::classify_web_environment(&[2.0, 1.5, 0.5], 0.0).unwrap();
    assert_eq!(env, cosmic_web::WebEnvironment::Node);

    let delta = cosmic_web::density_contrast(10.0, 1.0).unwrap();
    assert!(delta > 0.0, "nodes are overdense");
}

#[test]
fn test_sersic_elliptical_vs_spiral() {
    // n=4 (elliptical) drops off faster in center than n=1 (exponential disk)
    let i_elliptical = morphology::sersic_profile(0.1, 1.0, 4.0).unwrap();
    let i_spiral = morphology::sersic_profile(0.1, 1.0, 1.0).unwrap();
    // At R < R_e, elliptical (n=4) has steeper central concentration
    assert!(i_elliptical > i_spiral, "elliptical brighter in center");
}
