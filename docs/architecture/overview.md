# Architecture

## Module Map

- `src/cosmology.rs` — `Cosmology` struct (Planck 2018, WMAP9 presets), `FilterFunction` enum, A_s normalization cache via `OnceLock`
- `src/error.rs` — `BrahmandaError` enum with `thiserror`, input/output validation helpers (`require_finite`, `ensure_finite`)
- `src/constants.rs` — Physical constants only (C, G, hbar, k_B, M_sun, Mpc, kpc) — cosmological parameters live in `Cosmology`
- `src/morphology.rs` — Galaxy classification (Hubble types), Sersic profile, Faber-Jackson, Tully-Fisher, mass-metallicity (Tremonti + z evolution), Madau-Dickinson SFR density, cumulative stellar mass density
- `src/halo.rs` — NFW dark matter halos (density, enclosed mass, circular velocity), virial radius, concentration (Dutton-Maccio 2014), Press-Schechter / Sheth-Tormen / Tinker08 mass functions, Mo-White / ST / Tinker10 halo bias, mass accretion history (Wechsler 2002, McBride 2009), subhalo abundance matching (SHAM)
- `src/cosmic_web.rs` — Tidal tensor eigenvalue classification, density contrast, power-law correlation function, HSW void density profile, void-in-void excursion set (Sheth & van de Weygaert 2004), Hessian morphology classification, filamentarity, Minkowski functionals
- `src/power_spectrum.rs` — Eisenstein-Hu transfer function (no-wiggle + BAO), linear P(k) normalized to sigma_8, exact sigma(R) via adaptive Simpson integration with filter functions (top-hat, Gaussian, sharp-k), growth factor & rate, dark energy EoS (CPL), distances (comoving, luminosity, angular diameter), comoving volume element, lookback time, age of universe, angular C_l (Limber), Halofit nonlinear P(k), distance modulus, Sachs-Wolfe (ordinary + ISW)
- `src/logging.rs` — tracing-subscriber init (feature-gated behind `logging` feature)

## Central Design

All cosmological computations take `&Cosmology` as context. The `Cosmology` struct holds primary parameters (h, Omega_m, Omega_b, Omega_k, Omega_r, n_s, sigma_8, w0, wa, T_CMB), derived quantities (H0_SI, rho_crit, Omega_Lambda), and a lazily-cached power spectrum normalization amplitude A_s.

Pure mathematical functions (NFW profile, Sersic profile, Faber-Jackson, classification functions) take only their physics inputs — no cosmology context needed.

## Dependencies

- **hisab** — numerical integration (`integral_adaptive_simpson`), interpolation, linear algebra
- **serde** — serialization (all public types derive Serialize/Deserialize)
- **thiserror** — error enum derivation
- **tracing** — structured logging on error paths

## Integration: hisab-mimamsa Scale 5

brahmanda provides galactic density fields, cosmic web classification, and structure formation data for hisab-mimamsa's `unified::scale_bridge::bridge_scale_5()` (galactic structure to civilizational personality fields via bhava).
