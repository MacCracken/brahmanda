# Contributing to brahmanda

## Getting Started

```bash
git clone https://github.com/MacCracken/brahmanda.git
cd brahmanda
cargo test --all-features
```

## Guidelines

- All physics functions must validate inputs and return `Result` for fallible operations.
- Use `crate::constants` for physical constants — never define them locally.
- Use `&Cosmology` for all cosmological parameters — never hardcode Omega_m, sigma_8, etc.
- Tracing:
  - `error!` for computational failures (integration divergence, bisection failure)
  - `warn!` for domain boundary violations (unphysical redshift, NaN input)
  - `#[instrument]` on significant public functions (`level = "trace"` for computations, `level = "debug"` for aggregators)
  - Skip `#[instrument]` on trivial one-liners, tight-loop helpers, and hot constructors
- Every public function needs a doc comment with the formula it implements.
- Tests must validate against known physical values, not arbitrary numbers.
- Zero `unsafe`. Zero `println!`. Zero clippy warnings.
- `#[non_exhaustive]` on all public enums.
- `#[must_use]` on all pure functions returning computed results.
- `#[inline]` on hot-path functions.
- Serde derives on all public types.

## Testing

```bash
cargo test --all-features       # all tests
cargo test --no-default-features # no_std check
cargo clippy --all-features     # lint
cargo bench                     # benchmarks
make check                      # full quality gate
```

## Test Organization

- `tests/adversarial.rs` — NaN, Inf, edge case rejection
- `tests/invariants.rs` — physical invariants (monotonicity, bounds, scaling)
- `tests/serde_roundtrip.rs` — serialization roundtrip for all public types
- `tests/integration.rs` — cross-module integration and Display tests

## License

By contributing you agree that your contributions will be licensed under GPL-3.0-only.
