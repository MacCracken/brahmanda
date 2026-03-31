# Changelog

## 1.0.0 — 2026-03-31

### Added

- `angular_diameter_distance` — d_A(z) = χ(z)/(1+z)
- `comoving_volume_element` — dV/(dz dΩ) for survey volume calculations
- `lookback_time` — time elapsed since redshift z (Gyr)
- `age_of_universe` — cosmic age at redshift z (Gyr)
- `growth_rate` — f(z) = d ln D / d ln a, via numerical differentiation
- `tinker08_dndlnm` — Tinker et al. (2008) calibrated halo mass function (Δ=200)
- `bias_tinker10` — Tinker et al. (2010) calibrated halo bias (Δ=200)
- `Display` implementations for `HubbleType`, `WebEnvironment`, `HessianMorphology`
- `PartialEq` derives for `GalaxyProperties` and `HaloProperties`
- GitHub Actions CI workflow (fmt, clippy, test, deny, MSRV check)
- Benchmark runner script (`scripts/bench-history.sh`)
- Comprehensive test coverage: 243 tests (unit, integration, doc)

### Changed

- **BREAKING**: Version bump from 0.1.0 to 1.0.0 (stable API)
- Constants doc comments now cite specific sources (Planck 2018 TT,TE,EE+lowE+lensing, CODATA 2018, IAU)
- `classify_hessian_morphology` eigenvalue sort uses `total_cmp()` instead of `partial_cmp().unwrap()` for NaN safety
- `classify_hessian_morphology` eigenvalue sort uses fixed-size array instead of `Vec` allocation
- `sigma_r` documentation clarifies power-law approximation limitations

### Fixed

- Unresolved rustdoc link `[Mpc]` in `distance_modulus` doc comment

## 0.1.0 — 2026-03-31

### Added

- `morphology` — Galaxy classification (Hubble types), Sersic profiles, Faber-Jackson, Tully-Fisher,
  mass-metallicity relation (Tremonti 2004 + z evolution), Madau-Dickinson SFR, stellar mass density
- `halo` — NFW dark matter halo profiles, virial radius, concentration-mass (Dutton & Maccio 2014),
  circular velocity, Press-Schechter & Sheth-Tormen mass functions, Mo & White / ST bias,
  mass accretion history (Wechsler, McBride), subhalo abundance matching (SHAM)
- `cosmic_web` — Tidal tensor classification (node/filament/sheet/void), density contrast,
  correlation function, HSW void profile, void excursion set (Sheth & van de Weygaert 2004),
  Hessian morphology, filamentarity, Minkowski functionals
- `power_spectrum` — Eisenstein-Hu transfer function (no-wiggle + BAO), growth factor (ΛCDM + w₀w_a),
  σ(R), dark energy EoS (CPL), comoving/luminosity distance, distance modulus, angular C_l (Limber),
  Halofit nonlinear P(k), Sachs-Wolfe (ordinary + ISW)
- `constants` — Planck 2018 cosmological parameters, physical constants (SI)
- `error` — `BrahmandaError` with input/output validation helpers
- Integration tests, criterion benchmarks, runnable example
