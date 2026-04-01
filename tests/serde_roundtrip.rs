//! Serde roundtrip tests for all public types.

use brahmanda::Cosmology;
use brahmanda::cosmic_web::{HessianMorphology, WebEnvironment};
use brahmanda::halo::HaloProperties;
use brahmanda::morphology::{GalaxyProperties, HubbleType};

fn cosmo() -> Cosmology {
    Cosmology::planck2018()
}

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

#[test]
fn cosmology_roundtrip() {
    let c = cosmo();
    let json = serde_json::to_string(&c).unwrap();
    let back: Cosmology = serde_json::from_str(&json).unwrap();
    assert_eq!(c, back);
}
