# Brahmanda — Claude Code Instructions

## Project Identity

**Brahmanda** (Sanskrit: ब्रह्माण्ड — the cosmic egg, the universe) — Large-scale structure physics for AGNOS

- **Type**: Flat library crate
- **License**: GPL-3.0-only
- **MSRV**: 1.89
- **Version**: SemVer 1.0.0

## Consumers

hisab-mimamsa (unified::scale_bridge Scale 5), kiran/joshua (galaxy simulation), agnosai (cosmological reasoning)

## Architecture

- `src/cosmology.rs` — `Cosmology` struct (Planck 2018, WMAP9 presets), `FilterFunction` enum, A_s normalization cache
- `src/error.rs` — `BrahmandaError` enum with `thiserror`, input/output validation helpers
- `src/constants.rs` — Physical constants only (C, G, ℏ, k_B, M_sun, Mpc, kpc) — cosmological params live in `Cosmology`
- `src/morphology.rs` — Galaxy classification, Sersic profile, Faber-Jackson, Tully-Fisher, mass-metallicity, Madau-Dickinson SFR, stellar mass density
- `src/halo.rs` — NFW dark matter halos, virial radius, concentration, circular velocity, Press-Schechter / Sheth-Tormen / Tinker08 mass functions, Mo-White / ST / Tinker10 bias, mass accretion history (Wechsler, McBride), SHAM
- `src/cosmic_web.rs` — Tidal tensor classification, density contrast, correlation function, HSW void profile, void excursion set, Hessian morphology, filamentarity, Minkowski functionals
- `src/power_spectrum.rs` — Transfer function (no-wiggle + BAO), linear P(k), exact σ(R) with filter functions, growth factor & rate, dark energy w(z), distances (comoving, luminosity, angular diameter), comoving volume, lookback time, age of universe, angular C_l, Halofit nonlinear P(k), distance modulus, Sachs-Wolfe (ordinary + ISW)
- `src/logging.rs` — tracing-subscriber init (feature-gated)

## Dependencies

- **hisab** — numerical integration (adaptive Simpson), interpolation, linear algebra
- **serde** — serialization
- **thiserror** — error handling
- **tracing** — structured logging

## Roadmap

Future work lives in `docs/development/roadmap.md` (P7–P10 + quality track).

### Integration: hisab-mimamsa Scale 5

brahmanda provides galactic density fields, cosmic web classification, and structure formation data for hisab-mimamsa's `unified::scale_bridge::bridge_scale_5()` (galactic structure → civilizational personality fields via bhava).
