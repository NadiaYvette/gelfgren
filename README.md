# Gelfgren - Piecewise Rational Interpolation Library

A high-performance library for piecewise rational interpolation and approximation, written in Rust with bindings to C, C++, Java, Fortran, Python, R, Ruby, Haskell, and Mercury.

Named after Jan Gelfgren's 1975 paper "Piecewise Rational Interpolation", this library implements algorithms for constructing rational approximants (P_n/Q_m) on mesh intervals, with applications to boundary value problems and function approximation.

## Features

- **Bernstein Polynomial Representation**: Numerically stable polynomial operations
- **Rational Functions**: P_n(x)/Q_m(x) with automatic simplification
- **Padé Approximants**: From power series or derivative data
- **Lagrange-Hermite Interpolation**: Using Traub's formulas with Bell polynomials
- **Piecewise Rational Approximation**: Gelfgren's method for mesh-based construction
- **Boundary Value Problems**: Linear constraint systems for ODE solutions
- **Multi-Language FFI**: Bindings for 9+ programming languages

## Project Structure

```
gelfgren/
├── gelfgren-core/          # Core Rust implementation
├── gelfgren-ffi/           # C FFI layer
├── bindings/               # Language-specific bindings
│   ├── c/                  # C examples
│   ├── cpp/                # C++ wrappers
│   ├── java/               # JNI bindings
│   ├── fortran/            # ISO_C_BINDING
│   ├── python/             # PyO3 bindings
│   ├── r/                  # extendr bindings
│   ├── ruby/               # Magnus bindings
│   ├── haskell/            # FFI bindings
│   └── mercury/            # foreign_proc bindings
├── examples/               # Example programs
└── docs/                   # Documentation
```

## Building

### Requirements

- Rust 1.70+ (install from https://rustup.rs)
- Cargo

### Build Core Library

```bash
cargo build --release
```

### Run Tests

```bash
cargo test --workspace
```

## Theoretical Foundation

This library implements algorithms from three seminal papers:

1. **Gelfgren (1975)**: "Piecewise Rational Interpolation"
   - Core piecewise rational interpolation method
   - Symmetric Padé approximants on subintervals
   - Error estimation and convergence theory

2. **Traub (1964)**: "On Lagrange-Hermite Interpolation"
   - Interpolation with function and derivative values
   - Bell polynomial formulation
   - Foundation for Padé construction with derivative data

3. **Farouki & Rajan (1987)**: "Algorithms for Polynomials in Bernstein Form"
   - Numerically stable Bernstein polynomial representation
   - Complete set of polynomial operations
   - Proof of superior conditioning vs. power basis

## Development Status

- [x] Phase 0: Research and planning
- [ ] Phase 1: Bernstein polynomial foundation
- [ ] Phase 2: Rational functions
- [ ] Phase 3: Padé approximants
- [ ] Phase 4: Lagrange-Hermite interpolation
- [ ] Phase 5: Piecewise rational construction
- [ ] Phase 6: Boundary value problem support
- [ ] Phase 7: FFI layer (C bindings)
- [ ] Phase 8: Multi-language bindings
- [ ] Phase 9: CI/CD and documentation

See [IMPLEMENTATION_PLAN.md](docs/IMPLEMENTATION_PLAN.md) for detailed roadmap.

## License

MIT OR Apache-2.0

## Contributing

Contributions are welcome! Please see CONTRIBUTING.md for guidelines.
