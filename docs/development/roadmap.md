# Brahmanda Roadmap

Future development priorities for brahmanda. See CHANGELOG.md for what shipped.

## P7 — Halo Infrastructure

- [ ] Mass definition system: M200c, M200m, M500c, Mvir with Bryan-Norman Delta_vir(z) + conversions
- [ ] Einasto profile: rho(r), M(<r), alpha-mass relation (Navarro+ 2004, Gao+ 2008)
- [ ] Additional c(M,z) models: Diemer-Joyce 2019, Ishiyama+ 2021, Duffy+ 2008, Bullock+ 2001
- [ ] NFW Fourier transform u(k|M) — analytic Si/Ci form (Cooray & Sheth 2002)
- [ ] Splashback radius fitting functions (Diemer & Kravtsov 2014, More+ 2015)
- [ ] Additional mass functions: Watson 2013, Despali 2016, Bocquet 2016, Comparat 2017

## P8 — Halo Model & Galaxy Connection

- [ ] Halo model power spectrum: P_1h(k) + P_2h(k) (Cooray & Sheth 2002)
- [ ] HOD framework: Zheng+ 2005 5-parameter model (N_cen, N_sat)
- [ ] xi(r) from P(k) via numerical Hankel transform
- [ ] Redshift-space distortions: Kaiser P_s(k,mu) and multipoles P_0, P_2, P_4

## P9 — Observables & Survey Science

- [ ] Weak lensing convergence power spectrum C_l^kk with lensing efficiency kernel
- [ ] HMCode-2020 nonlinear P(k) (Mead+ 2021) — replaces Smith+ 2003 Halofit
- [ ] Curvature support (Omega_k != 0) in all distance functions
- [ ] Massive neutrino effects on P(k) and H(z) (N_eff, sum m_nu)
- [ ] Alcock-Paczynski geometric distortions

## P10 — Advanced / Niche

- [ ] Perturbation theory: SPT 1-loop P(k), EFTofLSS counter-terms
- [ ] Fisher matrix forecasting utilities
- [ ] Thermal SZ pressure profile (Arnaud+ 2010) and tSZ power spectrum
- [ ] Conditional mass function / conditional luminosity function
- [ ] Halo exclusion corrections for 2-halo term
- [ ] Excursion set merger rates and progenitor mass functions
- [ ] Lyman-alpha forest transmission F(z), tau_eff(z)

## Quality & Infrastructure

- [ ] Interpolation tables / caching for expensive integrals (growth factor, distance, sigma(R))
- [ ] Vectorized / batch interfaces for key functions (accept &[f64] -> Vec<f64>)
- [ ] Cross-validation tests against CCL/Colossus reference values at fixed tolerance
- [ ] Property-based testing with proptest (physical invariants across random inputs)
- [ ] no_std audit — verify no accidental std dependencies
