//! Basic brahmanda usage — halos, galaxies, cosmic web.

use brahmanda::{cosmic_web, halo, morphology, power_spectrum};

fn main() {
    // Dark matter halo for a Milky Way-mass galaxy
    let mw = halo::HaloProperties::from_mass(1e12).unwrap();
    println!("Milky Way dark matter halo:");
    println!("  Virial mass:   1e12 M_sun");
    println!("  Virial radius: {:.0} kpc", mw.r_vir_kpc);
    println!("  Scale radius:  {:.0} kpc", mw.r_s_kpc);
    println!("  Concentration: {:.1}", mw.concentration);

    // Galaxy morphology
    let sigma_mw = 100.0; // km/s velocity dispersion
    let sigma_giant = 300.0;
    let lum_ratio = morphology::faber_jackson_ratio(sigma_giant, sigma_mw).unwrap();
    println!("\nFaber-Jackson: σ={sigma_giant} vs σ={sigma_mw} → L ratio = {lum_ratio:.0}×");

    // Cosmic web classification
    let node = cosmic_web::classify_web_environment(&[2.0, 1.0, 0.5], 0.0).unwrap();
    let void = cosmic_web::classify_web_environment(&[-0.5, -0.8, -1.0], 0.0).unwrap();
    println!("\nCosmic web: [2,1,0.5] → {node:?}, [-0.5,-0.8,-1] → {void:?}");

    // Structure growth
    println!("\nGrowth factor D(z)/D(0):");
    for z in [0.0, 0.5, 1.0, 2.0, 5.0] {
        let d = power_spectrum::growth_factor(z, 0.315).unwrap();
        println!("  z={z:.1}: D = {d:.4}");
    }
}
