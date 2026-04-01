//! Physical invariant tests — monotonicity, bounds, scaling relations.

use brahmanda::Cosmology;
use brahmanda::cosmic_web;
use brahmanda::halo;
use brahmanda::morphology;
use brahmanda::power_spectrum;

fn cosmo() -> Cosmology {
    Cosmology::planck2018()
}

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
            "NFW density not decreasing: rho({})={} >= rho({})={}",
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
    assert!(s_small > s_8, "sigma(2) > sigma(8)");
    assert!(s_8 > s_large, "sigma(8) > sigma(32)");
}

#[test]
fn sigma_r_decreases_with_redshift() {
    let c = cosmo();
    let radii = [2.0, 8.0, 32.0];
    for &r in &radii {
        let s0 = power_spectrum::sigma_r(r, 0.0, &c).unwrap();
        let s1 = power_spectrum::sigma_r(r, 1.0, &c).unwrap();
        let s3 = power_spectrum::sigma_r(r, 3.0, &c).unwrap();
        assert!(s0 > s1, "sigma({r}, z=0) > sigma({r}, z=1)");
        assert!(s1 > s3, "sigma({r}, z=1) > sigma({r}, z=3)");
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
        assert!((xi - 1.0).abs() < 1e-10, "xi(r0) for gamma={gamma}: {xi}");
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
        "chi(z=1) = {chi}, expected ~3364 Mpc"
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
        "rho_SFR(z=0) = {sfr}, expected ~0.015"
    );
}
