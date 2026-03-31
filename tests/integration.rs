//! Cross-module integration tests.

use brahmanda::Cosmology;
use brahmanda::cosmic_web::{self, WebEnvironment};
use brahmanda::halo::{self, HaloProperties};
use brahmanda::morphology::{self, GalaxyProperties, HubbleType};
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

// ============================================================
// P1: Adversarial input tests — NaN, Inf, edge cases
// ============================================================

mod adversarial {
    use super::*;

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
}

// ============================================================
// P1: Serde roundtrip tests
// ============================================================

mod serde_roundtrip {
    use super::*;

    #[test]
    fn hubble_type_roundtrip() {
        let types = [
            HubbleType::Elliptical,
            HubbleType::Lenticular,
            HubbleType::Spiral,
            HubbleType::BarredSpiral,
            HubbleType::Irregular,
        ];
        for &t in &types {
            let json = serde_json::to_string(&t).unwrap();
            let back: HubbleType = serde_json::from_str(&json).unwrap();
            assert_eq!(t, back, "roundtrip failed for {t:?}");
        }
    }

    #[test]
    fn galaxy_properties_roundtrip() {
        let gp = GalaxyProperties {
            hubble_type: HubbleType::Spiral,
            stellar_mass_msun: 5e10,
            effective_radius_kpc: 3.5,
            sersic_index: 1.0,
            velocity_dispersion_km_s: 100.0,
        };
        let json = serde_json::to_string(&gp).unwrap();
        let back: GalaxyProperties = serde_json::from_str(&json).unwrap();
        assert_eq!(back.hubble_type, gp.hubble_type);
        assert!((back.stellar_mass_msun - gp.stellar_mass_msun).abs() < 1.0);
        assert!((back.effective_radius_kpc - gp.effective_radius_kpc).abs() < 1e-10);
        assert!((back.sersic_index - gp.sersic_index).abs() < 1e-10);
        assert!((back.velocity_dispersion_km_s - gp.velocity_dispersion_km_s).abs() < 1e-10);
    }

    #[test]
    fn web_environment_roundtrip() {
        let envs = [
            WebEnvironment::Node,
            WebEnvironment::Filament,
            WebEnvironment::Sheet,
            WebEnvironment::Void,
        ];
        for &e in &envs {
            let json = serde_json::to_string(&e).unwrap();
            let back: WebEnvironment = serde_json::from_str(&json).unwrap();
            assert_eq!(e, back, "roundtrip failed for {e:?}");
        }
    }

    #[test]
    fn halo_properties_roundtrip() {
        let c = cosmo();
        let halo = HaloProperties::from_mass(1e12, &c).unwrap();
        let json = serde_json::to_string(&halo).unwrap();
        let back: HaloProperties = serde_json::from_str(&json).unwrap();
        assert!((back.m_vir_msun - halo.m_vir_msun).abs() < 1.0);
        assert!((back.r_vir_kpc - halo.r_vir_kpc).abs() < 1e-6);
        assert!((back.r_s_kpc - halo.r_s_kpc).abs() < 1e-6);
        assert!((back.concentration - halo.concentration).abs() < 1e-10);
    }

    #[test]
    fn cosmology_roundtrip() {
        let c = cosmo();
        let json = serde_json::to_string(&c).unwrap();
        let back: Cosmology = serde_json::from_str(&json).unwrap();
        assert_eq!(c, back);
    }
}

// ============================================================
// P1: Physical invariant tests
// ============================================================

mod physical_invariants {
    use super::*;

    #[test]
    fn nfw_enclosed_mass_monotonic() {
        let rho_s = 1e7;
        let r_s = 20.0;
        let radii = [1.0, 5.0, 10.0, 50.0, 100.0, 500.0];
        let masses: Vec<f64> = radii
            .iter()
            .map(|&r| halo::nfw_enclosed_mass(r, rho_s, r_s).unwrap())
            .collect();
        for i in 1..masses.len() {
            assert!(
                masses[i] > masses[i - 1],
                "enclosed mass not monotonic: M({})={} <= M({})={}",
                radii[i],
                masses[i],
                radii[i - 1],
                masses[i - 1]
            );
        }
    }

    #[test]
    fn nfw_density_monotonic_decreasing() {
        let rho_s = 1e7;
        let r_s = 20.0;
        let radii = [1.0, 5.0, 10.0, 50.0, 100.0, 500.0];
        let densities: Vec<f64> = radii
            .iter()
            .map(|&r| halo::nfw_density(r, rho_s, r_s).unwrap())
            .collect();
        for i in 1..densities.len() {
            assert!(
                densities[i] < densities[i - 1],
                "NFW density not decreasing: ρ({})={} >= ρ({})={}",
                radii[i],
                densities[i],
                radii[i - 1],
                densities[i - 1]
            );
        }
    }

    #[test]
    fn growth_factor_monotonically_decreasing_with_z() {
        let c = cosmo();
        let redshifts = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0];
        let growths: Vec<f64> = redshifts
            .iter()
            .map(|&z| power_spectrum::growth_factor(z, &c).unwrap())
            .collect();
        for i in 1..growths.len() {
            assert!(
                growths[i] < growths[i - 1],
                "growth factor not decreasing: D(z={})={} >= D(z={})={}",
                redshifts[i],
                growths[i],
                redshifts[i - 1],
                growths[i - 1]
            );
        }
    }

    #[test]
    fn growth_factor_bounded() {
        let c = cosmo();
        let d0 = power_spectrum::growth_factor(0.0, &c).unwrap();
        assert!((d0 - 1.0).abs() < 1e-10);

        for z in [0.0, 1.0, 10.0, 100.0] {
            let d = power_spectrum::growth_factor(z, &c).unwrap();
            assert!(d > 0.0, "growth factor negative at z={z}: {d}");
        }
    }

    #[test]
    fn sigma_r_scales_with_radius() {
        let c = cosmo();
        let s_small = power_spectrum::sigma_r(2.0, 0.0, &c).unwrap();
        let s_8 = power_spectrum::sigma_r(8.0, 0.0, &c).unwrap();
        let s_large = power_spectrum::sigma_r(32.0, 0.0, &c).unwrap();
        assert!(s_small > s_8, "σ(2) > σ(8)");
        assert!(s_8 > s_large, "σ(8) > σ(32)");
    }

    #[test]
    fn sigma_r_decreases_with_redshift() {
        let c = cosmo();
        let radii = [2.0, 8.0, 32.0];
        for &r in &radii {
            let s0 = power_spectrum::sigma_r(r, 0.0, &c).unwrap();
            let s1 = power_spectrum::sigma_r(r, 1.0, &c).unwrap();
            let s3 = power_spectrum::sigma_r(r, 3.0, &c).unwrap();
            assert!(s0 > s1, "σ({r}, z=0) > σ({r}, z=1)");
            assert!(s1 > s3, "σ({r}, z=1) > σ({r}, z=3)");
        }
    }

    #[test]
    fn virial_radius_increases_with_mass() {
        let c = cosmo();
        let masses = [1e8, 1e10, 1e12, 1e14, 1e15];
        let radii: Vec<f64> = masses
            .iter()
            .map(|&m| halo::virial_radius(m, &c).unwrap())
            .collect();
        for i in 1..radii.len() {
            assert!(
                radii[i] > radii[i - 1],
                "virial radius not increasing: R({})={} <= R({})={}",
                masses[i],
                radii[i],
                masses[i - 1],
                radii[i - 1]
            );
        }
    }

    #[test]
    fn concentration_decreases_with_mass() {
        let c = cosmo();
        let masses = [1e10, 1e12, 1e14, 1e15];
        let concs: Vec<f64> = masses
            .iter()
            .map(|&m| halo::concentration_dutton_maccio(m, &c).unwrap())
            .collect();
        for i in 1..concs.len() {
            assert!(
                concs[i] < concs[i - 1],
                "concentration not decreasing: c({})={} >= c({})={}",
                masses[i],
                concs[i],
                masses[i - 1],
                concs[i - 1]
            );
        }
    }

    #[test]
    fn sersic_profile_at_r_e_is_unity() {
        for n in [0.5, 1.0, 2.0, 4.0, 8.0] {
            let i = morphology::sersic_profile(10.0, 10.0, n).unwrap();
            assert!(
                (i - 1.0).abs() < 0.05,
                "Sersic at R_e for n={n}: got {i}, expected ~1.0"
            );
        }
    }

    #[test]
    fn transfer_function_bounded() {
        let ks = [0.001, 0.01, 0.1, 1.0, 10.0];
        for &k in &ks {
            let t = power_spectrum::transfer_function_eh(k, 0.315, 0.049, 0.674).unwrap();
            assert!(t > 0.0 && t <= 1.01, "T(k={k}) = {t}, expected (0, 1]");
        }
    }

    #[test]
    fn correlation_function_at_r0_is_unity() {
        for gamma in [1.0, 1.5, 1.8, 2.0, 3.0] {
            let xi = cosmic_web::two_point_correlation_power_law(5.0, 5.0, gamma).unwrap();
            assert!((xi - 1.0).abs() < 1e-10, "ξ(r₀) for γ={gamma}: {xi}");
        }
    }

    #[test]
    fn void_radius_scales_as_cube_root() {
        let r1 = cosmic_web::void_radius_mpc(1000.0).unwrap();
        let r2 = cosmic_web::void_radius_mpc(2000.0).unwrap();
        let ratio = r2 / r1;
        let expected = 2.0_f64.cbrt();
        assert!(
            (ratio - expected).abs() < 1e-6,
            "void radius scaling: {ratio}, expected {expected}"
        );
    }

    #[test]
    fn comoving_distance_planck_z1() {
        let c = cosmo();
        let chi = power_spectrum::comoving_distance(1.0, &c).unwrap();
        assert!(
            (chi - 3364.0).abs() / 3364.0 < 0.05,
            "χ(z=1) = {chi}, expected ~3364 Mpc"
        );
    }

    #[test]
    fn luminosity_distance_z1() {
        let c = cosmo();
        let dl = power_spectrum::luminosity_distance(1.0, &c).unwrap();
        assert!(
            (dl - 6728.0).abs() / 6728.0 < 0.05,
            "d_L(z=1) = {dl}, expected ~6728 Mpc"
        );
    }

    #[test]
    fn nfw_circular_velocity_mw() {
        let c = cosmo();
        let halo = halo::HaloProperties::from_mass(1e12, &c).unwrap();
        let conc = halo.concentration;
        let r_s = halo.r_s_kpc;
        let g_c = (1.0 + conc).ln() - conc / (1.0 + conc);
        let rho_s = halo.m_vir_msun / (4.0 * std::f64::consts::PI * r_s.powi(3) * g_c);
        let v = halo::nfw_circular_velocity(8.0, rho_s, r_s).unwrap();
        assert!(
            v > 100.0 && v < 400.0,
            "MW v_c at 8 kpc: {v} km/s, expected ~220"
        );
    }

    #[test]
    fn press_schechter_value_range() {
        let c = cosmo();
        let n = halo::press_schechter_dndlnm(1e12, 0.0, &c).unwrap();
        assert!(n > 1e-6 && n < 1.0, "PS dn/dlnM(10^12, z=0) = {n}");
    }

    #[test]
    fn sfr_density_z0_value() {
        let sfr = morphology::sfr_density_madau_dickinson(0.0).unwrap();
        assert!(
            (sfr - 0.015).abs() / 0.015 < 0.3,
            "ρ_SFR(z=0) = {sfr}, expected ~0.015"
        );
    }
}

// ============================================================
// Display trait tests
// ============================================================

mod display_traits {
    use super::*;
    use brahmanda::cosmic_web::HessianMorphology;

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
}

// ============================================================
// Additional serde roundtrip
// ============================================================

mod serde_roundtrip_extended {
    use brahmanda::cosmic_web::HessianMorphology;

    #[test]
    fn hessian_morphology_roundtrip() {
        let variants = [
            HessianMorphology::Blob,
            HessianMorphology::Filament,
            HessianMorphology::Wall,
            HessianMorphology::Background,
        ];
        for &v in &variants {
            let json = serde_json::to_string(&v).unwrap();
            let back: HessianMorphology = serde_json::from_str(&json).unwrap();
            assert_eq!(v, back, "roundtrip failed for {v:?}");
        }
    }
}
