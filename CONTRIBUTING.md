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
- Add `tracing::warn!` on all physics-boundary error paths.
- Every public function needs a doc comment with the formula it implements.
- Tests must validate against known physical values, not arbitrary numbers.
- Zero `unsafe`. Zero `println!`. Zero clippy warnings.
- `#[non_exhaustive]` on all public enums.
- `#[must_use]` on all pure functions.
- `#[inline]` on hot-path functions.
- Serde derives on all public types.

## Testing

```bash
cargo test --all-features       # all tests
cargo clippy --all-features     # lint
cargo bench                     # benchmarks
```

## License

By contributing you agree that your contributions will be licensed under GPL-3.0-only.
