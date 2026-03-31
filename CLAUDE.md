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
- `src/halo.rs` — NFW dark matter halos, virial radius, concentration, circular velocity
- `src/cosmic_web.rs` — Tidal tensor classification, density contrast, correlation function
- `src/power_spectrum.rs` — Transfer function, growth factor, primordial power, σ(R)
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

- [ ] Adversarial input tests (NaN, Inf, edge cases)
- [ ] Serde roundtrip tests for all public types
- [ ] Physical invariant tests (mass conservation, growth monotonicity)
- [ ] Doc tests on all public functions

### P2 — Extensions

- [ ] Press-Schechter halo mass function
- [ ] Halo bias (Mo & White 1996)
- [ ] BAO wiggles in transfer function
- [ ] Void profile models (Hamaus et al.)
- [ ] Angular power spectrum C_l
- [ ] Dark energy equation of state w(z) models

### Integration: hisab-mimamsa Scale 5

brahmanda provides galactic density fields, cosmic web classification, and structure formation data for hisab-mimamsa's `unified::scale_bridge::bridge_scale_5()` (galactic structure → civilizational personality fields via bhava).
