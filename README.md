# Gelfgren - Multi-Language Numerical Computing Library

A high-performance numerical computing library written in Rust with bindings to C, C++, Java, Fortran, Python, R, Ruby, Haskell, and Mercury.

## Features

- **Linear Algebra**: Matrix and vector operations, decompositions, linear solvers
- **Statistics**: Descriptive statistics, distributions, hypothesis testing
- **Numerical Methods**: Integration, differentiation, ODE solvers, root finding, optimization
- **Special Functions**: Gamma, beta, Bessel functions

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

## Development Status

- [x] Phase 1: Rust core library with basic numerical primitives
- [ ] Phase 2: FFI layer with C bindings
- [ ] Phase 3: C/C++ bindings
- [ ] Phase 4: Python, Java, R bindings
- [ ] Phase 5: Ruby, Fortran, Haskell, Mercury bindings
- [ ] Phase 6: CI/CD pipeline
- [ ] Phase 7: Documentation
- [ ] Phase 8: Release automation
- [ ] Phase 9: Build orchestration scripts

## License

MIT OR Apache-2.0

## Contributing

Contributions are welcome! Please see CONTRIBUTING.md for guidelines.
