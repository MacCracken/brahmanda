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

- `src/error.rs` — `BrahmandaError` enum with `thiserror`, input/output validation helpers
- `src/constants.rs` — Planck 2018 cosmological parameters, physical constants (SI)
- `src/morphology.rs` — Galaxy classification, Sersic profile, Faber-Jackson, Tully-Fisher, mass-metallicity, Madau-Dickinson SFR, stellar mass density
- `src/halo.rs` — NFW dark matter halos, virial radius, concentration, circular velocity, Press-Schechter / Sheth-Tormen / Tinker08 mass functions, Mo-White / ST / Tinker10 bias, mass accretion history (Wechsler, McBride), SHAM
- `src/cosmic_web.rs` — Tidal tensor classification, density contrast, correlation function, HSW void profile, void excursion set, Hessian morphology, filamentarity, Minkowski functionals
- `src/power_spectrum.rs` — Transfer function (no-wiggle + BAO), growth factor & rate, primordial power, σ(R), dark energy w(z), distances (comoving, luminosity, angular diameter), comoving volume, lookback time, age of universe, angular C_l, Halofit nonlinear P(k), distance modulus, Sachs-Wolfe (ordinary + ISW)
- `src/logging.rs` — tracing-subscriber init (feature-gated)

## Dependencies

- **hisab** — trigonometry, linear algebra, numerical methods
- **serde** — serialization
- **thiserror** — error handling
- **tracing** — structured logging

## Roadmap

### P0 — Current (v0.1.0)

- [x] Galaxy morphology (Hubble types, Sersic, scaling relations)
- [x] Dark matter halos (NFW, virial, concentration)
- [x] Cosmic web topology (tidal tensor classification)
- [x] Power spectrum (transfer function, growth factor, σ(R))

### P1 — Hardening

- [x] Adversarial input tests (NaN, Inf, edge cases)
- [x] Serde roundtrip tests for all public types
- [x] Physical invariant tests (mass conservation, growth monotonicity)
- [x] Doc tests on all public functions

### P2 — Extensions

- [x] Press-Schechter halo mass function
- [x] Halo bias (Mo & White 1996)
- [x] BAO wiggles in transfer function
- [x] Void profile models (Hamaus et al.)
- [x] Angular power spectrum C_l
- [x] Dark energy equation of state w(z) models

### P3 — Halo Extensions

- [x] Sheth-Tormen mass function (ellipsoidal collapse)
- [x] Sheth-Tormen halo bias
- [x] Halo mass accretion history (Wechsler 2002, McBride 2009)
- [x] Subhalo abundance matching (Schechter SMF + SHAM)

### P4a — Power Spectrum / Cosmology Extensions

- [x] Nonlinear power spectrum (Smith et al. 2003 Halofit)
- [x] Luminosity distance and distance modulus μ(z)
- [x] Sachs-Wolfe effect (ordinary + integrated)

### P4b — Cosmic Web Extensions

- [x] Void-in-void excursion set (Sheth & van de Weygaert 2004)
- [x] Filament spine extraction (Hessian morphology, filamentarity)
- [x] Minkowski functionals for Gaussian excursion sets

### P4c — Galaxy Extensions

- [x] Mass-metallicity relation (Tremonti 2004 + redshift evolution)
- [x] Madau-Dickinson cosmic SFR density
- [x] Cumulative stellar mass density

### P5 — v1.0 Completeness

- [x] Angular diameter distance d_A(z)
- [x] Comoving volume element dV/dz
- [x] Lookback time and age of universe
- [x] Growth rate f(z) = d ln D / d ln a
- [x] Tinker 2008 mass function (calibrated, Δ=200)
- [x] Tinker 2010 halo bias (calibrated, Δ=200)

### Review Repairs (pre-publish)

- [x] Fix `cargo fmt` diff in `cosmic_web.rs:350` (array formatting)
- [x] Add doc comments to all constants in `constants.rs` — cite Planck 2018 / CODATA per-constant
- [x] Replace `partial_cmp().unwrap()` with `total_cmp()` in `cosmic_web.rs` eigenvalue sort for NaN safety
- [x] Add doc note on `sigma_r` power-law approximation — consumers expecting full integral should be aware
- [x] Display impls for all public enums
- [x] PartialEq for public structs with f64 fields
- [x] GitHub Actions CI workflow
- [x] Physics formula accuracy audit (all 10 core formulas verified)

### Integration: hisab-mimamsa Scale 5

brahmanda provides galactic density fields, cosmic web classification, and structure formation data for hisab-mimamsa's `unified::scale_bridge::bridge_scale_5()` (galactic structure → civilizational personality fields via bhava).
