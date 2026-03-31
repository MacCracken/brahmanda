# brahmanda

> **Brahmanda** (Sanskrit: ब्रह्माण्ड — the cosmic egg, the universe)

**Large-scale structure physics** — galaxy morphology, dark matter halos, cosmic web topology, and structure formation. Part of the [AGNOS](https://github.com/MacCracken) science stack, built on [hisab](https://crates.io/crates/hisab).

[![License: GPL-3.0-only](https://img.shields.io/badge/license-GPL--3.0--only-blue.svg)](LICENSE)

## Features

### Galaxy Morphology (`morphology`)

Hubble classification, Sersic surface brightness profiles, Faber-Jackson (elliptical) and Tully-Fisher (spiral) scaling relations.

### Dark Matter Halos (`halo`)

NFW density profiles, enclosed mass, circular velocity curves, virial radius, concentration-mass relations (Dutton & Maccio 2014).

### Cosmic Web (`cosmic_web`)

Tidal tensor eigenvalue classification (node/filament/sheet/void), density contrast, two-point correlation function, void radius estimation.

### Power Spectrum (`power_spectrum`)

Eisenstein-Hu transfer function, primordial power spectrum, linear growth factor (Carroll-Press-Turner 1992), σ(R) mass fluctuations.

## Quick Start

```rust
use brahmanda::halo::HaloProperties;

// Milky Way-mass dark matter halo
let halo = HaloProperties::from_mass(1e12).unwrap();
println!("Virial radius: {:.0} kpc", halo.r_vir_kpc);       // ~250 kpc
println!("Concentration: {:.1}", halo.concentration);         // ~8
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
