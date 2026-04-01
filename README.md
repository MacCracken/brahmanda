# brahmanda

> **Brahmanda** (Sanskrit: ब्रह्माण्ड — the cosmic egg, the universe)

**Large-scale structure physics** — galaxy morphology, dark matter halos, cosmic web topology, and structure formation. Part of the [AGNOS](https://github.com/MacCracken) science stack, built on [hisab](https://crates.io/crates/hisab).

[![License: GPL-3.0-only](https://img.shields.io/badge/license-GPL--3.0--only-blue.svg)](LICENSE)

## Features

### Cosmology (`cosmology`)

Parameterized cosmology with presets (Planck 2018, WMAP9). All computations take `&Cosmology` for full flexibility.

### Galaxy Morphology (`morphology`)

Hubble classification, Sersic surface brightness profiles, Faber-Jackson (elliptical) and Tully-Fisher (spiral) scaling relations, mass-metallicity relation, Madau-Dickinson SFR density.

### Dark Matter Halos (`halo`)

NFW density profiles, enclosed mass, circular velocity curves, virial radius, concentration-mass (Dutton & Maccio 2014). Mass functions (Press-Schechter, Sheth-Tormen, Tinker 2008), halo bias (Mo-White, ST, Tinker 2010), mass accretion history, SHAM.

### Cosmic Web (`cosmic_web`)

Tidal tensor eigenvalue classification (node/filament/sheet/void), density contrast, correlation function, void profiles, Hessian morphology, Minkowski functionals.

### Power Spectrum (`power_spectrum`)

Eisenstein-Hu transfer function (with BAO), exact sigma(R) integration, linear P(k) normalized to sigma_8, Halofit nonlinear P(k), growth factor & rate, dark energy EoS, distances (comoving, luminosity, angular diameter), lookback time, age of universe, angular C_l, Sachs-Wolfe.

## Quick Start

```rust
use brahmanda::{Cosmology, halo::HaloProperties};

let cosmo = Cosmology::planck2018();

// Milky Way-mass dark matter halo
let halo = HaloProperties::from_mass(1e12, &cosmo).unwrap();
println!("Virial radius: {:.0} kpc", halo.r_vir_kpc);       // ~250 kpc
println!("Concentration: {:.1}", halo.concentration);         // ~8
```

```rust
use brahmanda::{Cosmology, power_spectrum};

let cosmo = Cosmology::planck2018();

// Growth factor and age of the universe
let d = power_spectrum::growth_factor(1.0, &cosmo).unwrap();
let age = power_spectrum::age_of_universe(0.0, &cosmo).unwrap();
println!("D(z=1) = {d:.4}, age = {age:.1} Gyr");
```

```rust
use brahmanda::cosmic_web;

// Classify a cosmic web environment from tidal tensor eigenvalues
let env = cosmic_web::classify_web_environment(&[2.0, 1.0, 0.5], 0.0).unwrap();
assert_eq!(env, cosmic_web::WebEnvironment::Node);
```

## Relationship to AGNOS Science Stack

```
hisab (math foundation)
  ├── hisab-mimamsa — theoretical physics (GR, QFT, cosmology, unified)
  │     └── unified::scale_bridge — Scale 5 consumes brahmanda outputs
  ├── brahmanda (this) — galactic / large-scale structure
  ├── tara — stellar astrophysics
  ├── jyotish — astronomical computation
  └── falak — orbital mechanics
```

## License

GPL-3.0-only

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.
