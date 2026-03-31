# Brahmanda — Claude Code Instructions

## Project Identity

**Brahmanda** (Sanskrit: ब्रह्माण्ड — the cosmic egg, the universe) — Large-scale structure physics for AGNOS

- **Type**: Flat library crate
- **License**: GPL-3.0-only
- **MSRV**: 1.89
- **Version**: SemVer 0.1.0

## Consumers

hisab-mimamsa (unified::scale_bridge Scale 5), kiran/joshua (galaxy simulation), agnosai (cosmological reasoning)

## Architecture

- `src/error.rs` — `BrahmandaError` enum with `thiserror`, input/output validation helpers
- `src/constants.rs` — Planck 2018 cosmological parameters, physical constants (SI)
- `src/morphology.rs` — Galaxy classification, Sersic profile, Faber-Jackson, Tully-Fisher
- `src/halo.rs` — NFW dark matter halos, virial radius, concentration, circular velocity, Press-Schechter mass function, Mo & White bias
- `src/cosmic_web.rs` — Tidal tensor classification, density contrast, correlation function, HSW void profile
- `src/power_spectrum.rs` — Transfer function (no-wiggle + BAO), growth factor, primordial power, σ(R), dark energy w(z), comoving distance, angular power spectrum C_l
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

### Integration: hisab-mimamsa Scale 5

brahmanda provides galactic density fields, cosmic web classification, and structure formation data for hisab-mimamsa's `unified::scale_bridge::bridge_scale_5()` (galactic structure → civilizational personality fields via bhava).
