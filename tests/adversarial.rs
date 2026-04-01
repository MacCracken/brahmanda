//! Adversarial input tests — NaN, Inf, edge cases.

use brahmanda::Cosmology;
use brahmanda::cosmic_web;
use brahmanda::halo::{self, HaloProperties};
use brahmanda::morphology;
use brahmanda::power_spectrum;

fn cosmo() -> Cosmology {
    Cosmology::planck2018()
}

const POISON: [f64; 3] = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY];

// -- morphology --

#[test]
fn sersic_nan_inf() {
    for &bad in &POISON {
        assert!(morphology::sersic_profile(bad, 1.0, 4.0).is_err());
        assert!(morphology::sersic_profile(1.0, bad, 4.0).is_err());
        assert!(morphology::sersic_profile(1.0, 1.0, bad).is_err());
    }
}

#[test]
fn sersic_zero_and_negative() {
    assert!(morphology::sersic_profile(1.0, 0.0, 4.0).is_err());
    assert!(morphology::sersic_profile(1.0, -1.0, 4.0).is_err());
    assert!(morphology::sersic_profile(1.0, 1.0, 0.0).is_err());
    assert!(morphology::sersic_profile(1.0, 1.0, -1.0).is_err());
}

#[test]
fn sersic_r_zero_ok() {
    let result = morphology::sersic_profile(0.0, 1.0, 4.0);
    let _ = result;
}

#[test]
fn faber_jackson_nan_inf() {
    for &bad in &POISON {
        assert!(morphology::faber_jackson_ratio(bad, 200.0).is_err());
        assert!(morphology::faber_jackson_ratio(200.0, bad).is_err());
    }
}

#[test]
fn faber_jackson_zero_negative() {
    assert!(morphology::faber_jackson_ratio(0.0, 200.0).is_err());
    assert!(morphology::faber_jackson_ratio(200.0, 0.0).is_err());
    assert!(morphology::faber_jackson_ratio(-1.0, 200.0).is_err());
}

#[test]
fn tully_fisher_nan_inf() {
    for &bad in &POISON {
        assert!(morphology::tully_fisher_ratio(bad, 150.0, 4.0).is_err());
        assert!(morphology::tully_fisher_ratio(150.0, bad, 4.0).is_err());
        assert!(morphology::tully_fisher_ratio(150.0, 150.0, bad).is_err());
    }
}

// -- halo --

#[test]
fn nfw_density_nan_inf() {
    for &bad in &POISON {
        assert!(halo::nfw_density(bad, 1e7, 20.0).is_err());
        assert!(halo::nfw_density(10.0, bad, 20.0).is_err());
        assert!(halo::nfw_density(10.0, 1e7, bad).is_err());
    }
}

#[test]
fn nfw_density_zero_negative() {
    assert!(halo::nfw_density(0.0, 1e7, 20.0).is_err());
    assert!(halo::nfw_density(10.0, 0.0, 20.0).is_err());
    assert!(halo::nfw_density(10.0, 1e7, 0.0).is_err());
    assert!(halo::nfw_density(-5.0, 1e7, 20.0).is_err());
}

#[test]
fn nfw_enclosed_mass_nan_inf() {
    for &bad in &POISON {
        assert!(halo::nfw_enclosed_mass(bad, 1e7, 20.0).is_err());
        assert!(halo::nfw_enclosed_mass(10.0, bad, 20.0).is_err());
        assert!(halo::nfw_enclosed_mass(10.0, 1e7, bad).is_err());
    }
}

#[test]
fn nfw_circular_velocity_nan_inf() {
    for &bad in &POISON {
        assert!(halo::nfw_circular_velocity(bad, 1e7, 20.0).is_err());
        assert!(halo::nfw_circular_velocity(10.0, bad, 20.0).is_err());
        assert!(halo::nfw_circular_velocity(10.0, 1e7, bad).is_err());
    }
}

#[test]
fn virial_radius_nan_inf() {
    let c = cosmo();
    for &bad in &POISON {
        assert!(halo::virial_radius(bad, &c).is_err());
    }
    assert!(halo::virial_radius(0.0, &c).is_err());
    assert!(halo::virial_radius(-1e12, &c).is_err());
}

#[test]
fn concentration_nan_inf() {
    let c = cosmo();
    for &bad in &POISON {
        assert!(halo::concentration_dutton_maccio(bad, &c).is_err());
    }
    assert!(halo::concentration_dutton_maccio(0.0, &c).is_err());
    assert!(halo::concentration_dutton_maccio(-1e12, &c).is_err());
}

#[test]
fn halo_properties_nan_inf() {
    let c = cosmo();
    for &bad in &POISON {
        assert!(HaloProperties::from_mass(bad, &c).is_err());
    }
}

// -- cosmic_web --

#[test]
fn classify_web_nan_inf() {
    for &bad in &POISON {
        assert!(cosmic_web::classify_web_environment(&[bad, 1.0, 0.5], 0.0).is_err());
        assert!(cosmic_web::classify_web_environment(&[1.0, 1.0, 0.5], bad).is_err());
    }
}

#[test]
fn density_contrast_nan_inf() {
    for &bad in &POISON {
        assert!(cosmic_web::density_contrast(bad, 1.0).is_err());
        assert!(cosmic_web::density_contrast(1.0, bad).is_err());
    }
    assert!(cosmic_web::density_contrast(1.0, 0.0).is_err());
    assert!(cosmic_web::density_contrast(1.0, -1.0).is_err());
}

#[test]
fn void_radius_nan_inf() {
    for &bad in &POISON {
        assert!(cosmic_web::void_radius_mpc(bad).is_err());
    }
    assert!(cosmic_web::void_radius_mpc(0.0).is_err());
    assert!(cosmic_web::void_radius_mpc(-100.0).is_err());
}

#[test]
fn two_point_correlation_nan_inf() {
    for &bad in &POISON {
        assert!(cosmic_web::two_point_correlation_power_law(bad, 5.0, 1.8).is_err());
        assert!(cosmic_web::two_point_correlation_power_law(5.0, bad, 1.8).is_err());
        assert!(cosmic_web::two_point_correlation_power_law(5.0, 5.0, bad).is_err());
    }
}

// -- power_spectrum --

#[test]
fn transfer_function_nan_inf() {
    for &bad in &POISON {
        assert!(power_spectrum::transfer_function_eh(bad, 0.315, 0.049, 0.674).is_err());
        assert!(power_spectrum::transfer_function_eh(0.1, bad, 0.049, 0.674).is_err());
        assert!(power_spectrum::transfer_function_eh(0.1, 0.315, bad, 0.674).is_err());
        assert!(power_spectrum::transfer_function_eh(0.1, 0.315, 0.049, bad).is_err());
    }
}

#[test]
fn transfer_function_zero_negative() {
    assert!(power_spectrum::transfer_function_eh(0.0, 0.315, 0.049, 0.674).is_err());
    assert!(power_spectrum::transfer_function_eh(-1.0, 0.315, 0.049, 0.674).is_err());
    assert!(power_spectrum::transfer_function_eh(0.1, 0.0, 0.049, 0.674).is_err());
    assert!(power_spectrum::transfer_function_eh(0.1, 0.315, 0.049, 0.0).is_err());
}

#[test]
fn primordial_power_nan_inf() {
    let c = cosmo();
    for &bad in &POISON {
        assert!(power_spectrum::primordial_power(bad, &c).is_err());
    }
    assert!(power_spectrum::primordial_power(0.0, &c).is_err());
    assert!(power_spectrum::primordial_power(-1.0, &c).is_err());
}

#[test]
fn growth_factor_nan_inf() {
    let c = cosmo();
    for &bad in &POISON {
        assert!(power_spectrum::growth_factor(bad, &c).is_err());
    }
    assert!(power_spectrum::growth_factor(-2.0, &c).is_err());
}

#[test]
fn sigma_r_nan_inf() {
    let c = cosmo();
    for &bad in &POISON {
        assert!(power_spectrum::sigma_r(bad, 0.0, &c).is_err());
    }
    assert!(power_spectrum::sigma_r(0.0, 0.0, &c).is_err());
    assert!(power_spectrum::sigma_r(-1.0, 0.0, &c).is_err());
}

// -- extreme but valid inputs --

#[test]
fn extreme_halo_masses() {
    let c = cosmo();
    let dwarf = HaloProperties::from_mass(1e8, &c).unwrap();
    assert!(dwarf.r_vir_kpc > 0.0);
    assert!(dwarf.concentration > 0.0);

    let cluster = HaloProperties::from_mass(1e15, &c).unwrap();
    assert!(cluster.r_vir_kpc > dwarf.r_vir_kpc);
    assert!(cluster.concentration > 0.0);
}

#[test]
fn extreme_sersic_indices() {
    let i_low = morphology::sersic_profile(1.0, 1.0, 0.5).unwrap();
    assert!(i_low.is_finite());

    let i_high = morphology::sersic_profile(1.0, 1.0, 10.0).unwrap();
    assert!(i_high.is_finite());
}

#[test]
fn growth_factor_boundary() {
    let c = cosmo();
    let d = power_spectrum::growth_factor(-1.0, &c);
    let _ = d;
}

#[test]
fn very_high_redshift() {
    let c = cosmo();
    let d = power_spectrum::growth_factor(100.0, &c).unwrap();
    assert!(d > 0.0 && d < 0.1, "growth factor at z=100: {d}");
}
