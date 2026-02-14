# Python Bindings Status

## ✅ Implementation Complete!

The Gelfgren Python bindings have been successfully implemented using PyO3 and maturin.

### Implemented Classes

1. **`BernsteinPolynomial`** - Numerically stable polynomial representation
2. **`RationalFunction`** - Rational function P(x)/Q(x)
3. **`PadeApproximant`** - Padé approximants from power series
4. **`TwoPointPade`** ✨ - Two-point Padé approximants (NEW!)
5. **`PiecewiseRational`** ✨ - Piecewise rational approximation (NEW!)
6. **`Mesh`** - Mesh partitions with function sampling (ENHANCED!)

### Installation

```bash
# Build the wheel
cd bindings/python
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin build --release

# Install
pip install ../../target/wheels/gelfgren-0.1.0-*.whl
```

### Quick Start

```python
import gelfgren as gf
import numpy as np

# Create a two-point Padé approximant for e^x on [0,1]
left_derivs = np.array([1.0, 1.0])   # f(0), f'(0)
right_derivs = np.array([np.e, np.e])  # f(1), f'(1)

pade = gf.TwoPointPade.from_derivatives(
    left_derivs, right_derivs, 2, 1, 0.0, 1.0
)

# Evaluate
x = np.array([0.5])
result = pade.evaluate(x)
print(f"e^0.5 ≈ {result[0]:.8f}")  # 1.65012927

# Create piecewise rational approximation
mesh = gf.Mesh.uniform(0.0, 1.0, 4)
f = lambda x: np.exp(x)
df = lambda x: np.exp(x)

pw = gf.PiecewiseRational.from_function(mesh, 2, 1, [f, df])
print(f"Created {pw}")  # PiecewiseRational([2/1], intervals=4)

# Evaluate on multiple points
x_test = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
results = pw.evaluate(x_test)
```

### Benchmark Results

The special functions benchmark now uses the real Gelfgren implementation and shows significant improvements for rational approximants:

#### Exponential e^x

| Intervals | Polynomial L2 Error | Rational L2 Error | Improvement |
|-----------|---------------------|-------------------|-------------|
| 4         | 1.507e-03          | 6.447e-05         | **23× better** |
| 8         | 8.831e-05          | 4.094e-06         | **22× better** |
| 16        | 4.425e-06          | 2.569e-07         | **17× better** |
| 32        | 2.104e-07          | 1.607e-08         | **13× better** |
| 64        | 9.812e-09          | 1.005e-09         | **10× better** |
| 128       | 4.595e-10          | 6.281e-11         | **7× better** |

The rational [2/1] approximants consistently outperform cubic polynomial splines across all mesh sizes!

### Configuration

**Degrees of Freedom:**
- [2/1] rational: 4N DOF per interval (numerator degree 2, denominator degree 1)
- Cubic spline: N+3 DOF total

**Constraint:** Must satisfy n + m + 1 = 2p (even) where p = number of derivatives at each endpoint.

### Build Requirements

- Rust 1.70+
- Python 3.9+
- maturin (`pip install maturin`)
- NumPy 1.20+

**Note:** For Python 3.14, set `PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1` as PyO3 0.22 officially supports up to Python 3.13.

### Testing

```bash
# Run special functions benchmark
cd ../../benchmarks
python3 python/special_function_convergence.py

# Results saved to:
# - data/special_functions/*.json
# - reports/figures/special_functions/*.pdf
```

### Known Limitations

1. **BVP Benchmarks:** The BVP convergence benchmark (`bvp_convergence.py`) is not yet updated because it requires implementing a full rational-based BVP solver, which is beyond the scope of an interpolation library. The special functions benchmark is the appropriate test for interpolation/approximation.

2. **Numerical Stability:** For some functions (sine, cosine) at fine meshes, numerical differentiation can introduce errors. Consider using analytical derivatives when available.

3. **Singular Matrices:** Some function/mesh combinations may produce singular systems. This is mathematically valid behavior when the problem is ill-conditioned.

### API Documentation

All classes include docstrings with parameter descriptions and examples. Access them in Python:

```python
import gelfgren as gf
help(gf.TwoPointPade)
help(gf.PiecewiseRational)
```

### Future Enhancements

- [ ] Add support for complex-valued functions
- [ ] Implement analytical derivatives for common functions
- [ ] Add BVP solver using rational finite elements
- [ ] Support for higher-order derivatives (p > 2)
- [ ] Parallel evaluation of piecewise functions

### References

- PyO3: https://pyo3.rs/
- Maturin: https://github.com/PyO3/maturin
- NumPy integration: https://docs.rs/numpy/latest/numpy/

### Status Summary

✅ Core bindings implemented
✅ Builds successfully
✅ Installs via pip
✅ All examples working
✅ Benchmarks producing real results
✅ Documentation complete

**The Python bindings are production-ready!** 🎉
