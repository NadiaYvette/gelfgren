# Gelfgren Python Bindings

Python bindings for the Gelfgren piecewise rational interpolation library.

## Installation

### From source (development)

```bash
# Install maturin
pip install maturin

# Build and install in development mode
maturin develop --release

# Or build a wheel
maturin build --release
pip install target/wheels/gelfgren-*.whl
```

## Quick Start

```python
import gelfgren
import numpy as np

# Create a Bernstein polynomial
p = gelfgren.BernsteinPolynomial([1.0, 2.0, 3.0], 0.0, 1.0)

# Evaluate at points
x = np.array([0.0, 0.5, 1.0])
y = p.evaluate(x)

# Compute derivative
p_prime = p.derivative()

# Polynomial arithmetic
q = gelfgren.BernsteinPolynomial([1.0, 1.0], 0.0, 1.0)
sum_poly = p + q
```

## Features

- **Bernstein Polynomials**: Numerically stable polynomial operations
- **Rational Functions**: P(x)/Q(x) with automatic pole detection
- **Padé Approximants**: Rational approximations from power series
- **NumPy Integration**: Vectorized evaluation over arrays
- **Exception Handling**: Pythonic error reporting

## Documentation

See `example.py` for comprehensive usage examples.

## Requirements

- Python ≥ 3.9
- NumPy ≥ 1.20
- Rust (for building from source)

## License

MIT OR Apache-2.0
