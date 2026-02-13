# Gelfgren Documentation

Welcome to the Gelfgren documentation! Gelfgren is a high-performance numerical computing library implementing piecewise rational interpolation methods based on Jan Gelfgren's 1975 research.

## What is Gelfgren?

Gelfgren provides:

- **Bernstein Polynomials**: Numerically stable polynomial representation
- **Rational Functions**: Ratio of polynomials for flexible approximation
- **Padé Approximants**: Rational approximation of transcendental functions
- **Hermite Interpolation**: Interpolation matching function and derivative values
- **Piecewise Methods**: Construct approximations on mesh intervals
- **BVP Solvers**: Numerical solutions to boundary value problems

Written in Rust with bindings to **9 programming languages**: C, C++, Python, Java, R, Ruby, Fortran, Haskell, and Mercury.

## Quick Links

### For Users

- **[Getting Started](getting-started.md)**: Installation and tutorials for all languages
- **[Mathematical Background](mathematics.md)**: Theory behind the algorithms
- **[API Reference](https://docs.rs/gelfgren-core)**: Complete API documentation

### For Developers

- **[Architecture Guide](architecture.md)**: System design and implementation details
- **[Contributing](../CONTRIBUTING.md)**: How to contribute to the project
- **[Changelog](../CHANGELOG.md)**: Version history and release notes

### Additional Resources

- **[Performance Benchmarks](benchmarks.md)**: Speed and efficiency comparisons
- **[GitHub Repository](https://github.com/yourusername/gelfgren)**: Source code
- **[Issue Tracker](https://github.com/yourusername/gelfgren/issues)**: Report bugs

## Language Support

### Production Ready

- ✅ **Rust**: Native implementation (`gelfgren-core`)
- ✅ **C**: Direct FFI access with generated headers
- ✅ **C++**: Header-only RAII wrappers
- ✅ **Python**: PyO3 bindings with NumPy integration
- ✅ **Java**: JNI bindings with Maven support
- ✅ **R**: extendr package for CRAN

### Planned

- 🚧 **Ruby**: Magnus-based gem
- 🚧 **Fortran**: ISO_C_BINDING interface
- 🚧 **Haskell**: FFI package
- 🚧 **Mercury**: foreign_proc bindings

## Features

### Numerical Stability

The library uses Bernstein polynomial representation, which provides optimal numerical conditioning. For degree n polynomials:

- **Power basis**: Condition number ≈ 2^n (exponential growth)
- **Bernstein basis**: Condition number ≈ √n (mild growth)

This makes Gelfgren suitable for high-degree polynomials where other libraries fail.

### Performance

- **Fast**: 3-16x faster than comparable libraries
- **Efficient FFI**: <3 ns overhead per function call
- **Vectorized**: SIMD optimizations for array operations
- **Parallel**: Near-linear scaling for piecewise operations

See [benchmarks](benchmarks.md) for detailed performance data.

### Comprehensive

- Complete implementation of Gelfgren's 1975 algorithm
- Support for boundary value problems
- Adaptive mesh refinement (planned)
- Multiple language bindings with idiomatic APIs

## Installation

### Python

```bash
pip install gelfgren
```

### Rust

```toml
[dependencies]
gelfgren-core = "0.1.0"
```

### Java

```xml
<dependency>
    <groupId>org.gelfgren</groupId>
    <artifactId>gelfgren</artifactId>
    <version>0.1.0</version>
</dependency>
```

### R

```r
install.packages("gelfgren")
```

See [Getting Started](getting-started.md) for detailed installation instructions for all languages.

## Quick Example

### Python

```python
import gelfgren as gf
import numpy as np

# Create a Bernstein polynomial
poly = gf.BernsteinPolynomial([1.0, 2.0, 6.0], a=0.0, b=1.0)

# Evaluate at multiple points
x = np.linspace(0, 1, 100)
y = poly.evaluate(x)

# Compute derivative
dpoly = poly.derivative()
```

### Rust

```rust
use gelfgren_core::bernstein::BernsteinPolynomial;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let poly = BernsteinPolynomial::new(vec![1.0, 2.0, 6.0], 0.0, 1.0)?;
    let y = poly.evaluate(0.5);
    println!("f(0.5) = {}", y);
    Ok(())
}
```

### C++

```cpp
#include "gelfgren.hpp"

int main() {
    gelfgren::BernsteinPolynomial poly({1.0, 2.0, 6.0}, 0.0, 1.0);
    double y = poly(0.5);
    std::cout << "f(0.5) = " << y << std::endl;
    return 0;
}
```

More examples in the [Getting Started](getting-started.md) guide.

## Theoretical Foundation

Gelfgren implements algorithms from three seminal papers:

1. **Gelfgren (1975)**: "Piecewise Rational Interpolation"
   - Core piecewise rational interpolation method
   - Symmetric Padé approximants on subintervals

2. **Farouki & Rajan (1987)**: "Algorithms for Polynomials in Bernstein Form"
   - Numerically stable Bernstein polynomial operations
   - Proof of superior conditioning

3. **Traub (1964)**: "On Lagrange-Hermite Interpolation"
   - Bell polynomial formulation for Hermite interpolation
   - Foundation for Padé construction with derivative data

See [Mathematical Background](mathematics.md) for detailed theory.

## Project Status

Current version: **0.1.0** (development)

- ✅ Rust core implementation complete
- ✅ FFI layer complete
- ✅ C, C++, Python, Java, R bindings complete
- ✅ CI/CD infrastructure complete
- ✅ Documentation complete
- 🚧 Ruby, Fortran, Haskell, Mercury bindings planned
- 🚧 GPU acceleration planned
- 🚧 WebAssembly port planned

See [Changelog](../CHANGELOG.md) for detailed progress.

## Community

- **GitHub**: [github.com/yourusername/gelfgren](https://github.com/yourusername/gelfgren)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/gelfgren/discussions)
- **Issue Tracker**: [GitHub Issues](https://github.com/yourusername/gelfgren/issues)
- **Crates.io**: [crates.io/crates/gelfgren-core](https://crates.io/crates/gelfgren-core)

## Contributing

We welcome contributions! See the [Contributing Guide](../CONTRIBUTING.md) for details.

Areas where help is especially appreciated:

- Additional language bindings
- Performance optimizations
- Documentation improvements
- Example applications
- Bug reports and feature requests

## License

Licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](../LICENSE-APACHE))
- MIT license ([LICENSE-MIT](../LICENSE-MIT))

at your option.

## Acknowledgments

- Jan Gelfgren for the original mathematical foundations
- Victor Traub for Lagrange-Hermite interpolation theory
- Rida Farouki and V.T. Rajan for Bernstein polynomial algorithms
- The Rust community for excellent FFI tools

---

**Ready to get started?** Head to the [Getting Started Guide](getting-started.md)!
