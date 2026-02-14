# Hermite Interpolation with Rational Approximants

## Overview

**Yes!** Gelfgren supports **Hermite-type interpolation** for rational approximants, where you can interpolate not just function values but also **first, second, and higher-order derivatives**.

This is the core capability of **two-point Padé approximants**.

## What is Hermite Interpolation?

### Classical Hermite Interpolation (Polynomial)

Given data at n+1 points:
- **Function values**: f(x₀), f(x₁), ..., f(xₙ)
- **First derivatives**: f'(x₀), f'(x₁), ..., f'(xₙ)
- **Higher derivatives** (optional): f''(xᵢ), f'''(xᵢ), etc.

Find a polynomial P(x) such that:
```
P(xᵢ) = f(xᵢ)
P'(xᵢ) = f'(xᵢ)
P''(xᵢ) = f''(xᵢ)  (if provided)
...
```

**Degrees of freedom**: If you have p derivatives at each of n+1 points, you need a polynomial of degree p(n+1) - 1.

### Rational Hermite Interpolation (Gelfgren)

Same idea, but with **rational functions** R(x) = P(x)/Q(x):
```
R(xᵢ) = f(xᵢ)
R'(xᵢ) = f'(xᵢ)
R''(xᵢ) = f''(xᵢ)
...
```

**Advantage**: Rationals can represent poles, asymptotes, and other features that polynomials cannot.

## Two-Point Padé: The Core Implementation

### API

```rust
TwoPointPade::from_endpoint_derivatives(
    left_derivatives: &[f64],   // [f(a), f'(a), f''(a), ...]
    right_derivatives: &[f64],  // [f(b), f'(b), f''(b), ...]
    n: usize,                   // Numerator degree
    m: usize,                   // Denominator degree
    a: f64,                     // Left endpoint
    b: f64,                     // Right endpoint
)
```

### How Many Derivatives?

The number of derivatives p is determined by the constraint:
```
n + m + 1 = 2p  (must be even!)
```

So:
- **[2/1] rational**: n=2, m=1 → 2+1+1=4 → **p=2** (need f, f' at each endpoint)
- **[3/2] rational**: n=3, m=2 → 3+2+1=6 → **p=3** (need f, f', f'' at each endpoint)
- **[4/3] rational**: n=4, m=3 → 4+3+1=8 → **p=4** (need f, f', f'', f''' at each endpoint)

### Example: Interpolate e^x with First Derivatives

```python
import gelfgren as gf
import numpy as np

# Given: f(0)=1, f'(0)=1, f(1)=e, f'(1)=e
left_derivs = np.array([1.0, 1.0])          # [f(0), f'(0)]
right_derivs = np.array([np.e, np.e])        # [f(1), f'(1)]

# Create [2/1] rational matching these 4 conditions
pade = gf.TwoPointPade.from_derivatives(
    left_derivs, right_derivs,
    n=2, m=1,              # [2/1] rational
    a=0.0, b=1.0
)

# The rational now satisfies:
# R(0) = 1, R'(0) = 1, R(1) = e, R'(1) = e
print(pade.eval_scalar(0.5))  # ≈ 1.65012927 (actual: 1.64872127)
```

### Example: Interpolate with Second Derivatives

```python
# Given: f, f', f'' at both endpoints
# For e^x: all derivatives equal the function value
left_derivs = np.array([1.0, 1.0, 1.0])      # [f(0), f'(0), f''(0)]
right_derivs = np.array([np.e, np.e, np.e])  # [f(1), f'(1), f''(1)]

# Create [3/2] rational matching these 6 conditions
pade = gf.TwoPointPade.from_derivatives(
    left_derivs, right_derivs,
    n=3, m=2,              # [3/2] rational (3+2+1=6, p=3)
    a=0.0, b=1.0
)

# Even more accurate!
print(pade.eval_scalar(0.5))  # ≈ 1.64872128 (almost exact)
```

## Degrees of Freedom Breakdown

### [2/1] Rational with p=2

**Constraints** (4 total):
1. R(a) = f(a)
2. R'(a) = f'(a)
3. R(b) = f(b)
4. R'(b) = f'(b)

**Unknowns**:
- Numerator P(x): degree 2 → 3 coefficients
- Denominator Q(x): degree 1 → 2 coefficients
- Normalization: b₀ + b₁ = 1 removes 1 DOF
- **Total**: 3 + 2 - 1 = **4 DOF** ✓

### [3/2] Rational with p=3

**Constraints** (6 total):
1. R(a) = f(a)
2. R'(a) = f'(a)
3. R''(a) = f''(a)
4. R(b) = f(b)
5. R'(b) = f'(b)
6. R''(b) = f''(b)

**Unknowns**:
- Numerator: degree 3 → 4 coefficients
- Denominator: degree 2 → 3 coefficients
- Normalization: removes 1 DOF
- **Total**: 4 + 3 - 1 = **6 DOF** ✓

## Mathematical Details

### Derivative Matching Conditions

For a rational R(x) = P(x)/Q(x), the derivatives are:

```
R'(x) = [P'(x)Q(x) - P(x)Q'(x)] / Q(x)²

R''(x) = [P''(x)Q(x)² - 2P'(x)Q(x)Q'(x) - P(x)Q'(x)² + 2P(x)Q(x)Q''(x)] / Q(x)³
```

To match R^(k)(xᵢ) = f^(k)(xᵢ), we set up equations involving P and Q coefficients.

### Linear System

The implementation converts these conditions into a linear system:

**For k-th derivative at left endpoint**:
```
Σⱼ₌₀^k C(k,j) · [P^(j)(a) · Q^(k-j)(a)] = f^(k)(a) · Q^(k)(a)
```

Using **forward difference formulas** for derivatives in Bernstein basis:
```
P^(k)(a) = k!/(n-k)! · 1/Δx^k · Σⱼ₌₀^k (-1)^(k-j) C(k,j) aⱼ
```

**For k-th derivative at right endpoint**:

Using **backward difference formulas**:
```
P^(k)(b) = k!/(n-k)! · 1/Δx^k · Σⱼ₌₀^k (-1)^j C(k,j) aₙ₋ⱼ
```

These formulas are implemented in `pade/two_point.rs` in the `construct_rational_pade()` function.

## Piecewise Hermite Interpolation

### Multiple Intervals

The `PiecewiseRational` class applies two-point Padé to multiple subintervals:

```python
import gelfgren as gf
import numpy as np

# Create mesh: [0, 0.25, 0.5, 0.75, 1.0]
mesh = gf.Mesh.uniform(0.0, 1.0, 4)

# Sample function and first derivative
f = lambda x: np.exp(x)
df = lambda x: np.exp(x)

# This creates 4 subintervals, each with [2/1] rational
# matching f and f' at subinterval endpoints
pw = gf.PiecewiseRational.from_function(mesh, 2, 1, [f, df])

# Evaluate anywhere in [0,1]
x_test = np.linspace(0, 1, 100)
y_approx = pw.evaluate(x_test)
```

### Automatic Derivative Computation

If you only have the function (not derivatives), you can compute derivatives numerically:

```python
def compute_derivatives(f, x, order):
    """Compute derivatives up to given order using finite differences."""
    h = 1e-7
    derivs = [f(x)]  # f(x)

    if order >= 1:
        # First derivative: [f(x+h) - f(x-h)] / 2h
        derivs.append((f(x + h) - f(x - h)) / (2*h))

    if order >= 2:
        # Second derivative: [f(x+h) - 2f(x) + f(x-h)] / h²
        derivs.append((f(x + h) - 2*f(x) + f(x - h)) / (h*h))

    return derivs

# Use for each mesh point
left_derivs = compute_derivatives(f, 0.0, order=1)
right_derivs = compute_derivatives(f, 1.0, order=1)

pade = gf.TwoPointPade.from_derivatives(
    left_derivs, right_derivs, 2, 1, 0.0, 1.0
)
```

## Comparison: Polynomial vs Rational Hermite

### Polynomial Hermite (Cubic)

**With f and f' at 2 points** (4 conditions):
- Need degree 3 polynomial
- 4 DOF from 4 coefficients
- Good for smooth functions
- Cannot represent poles or asymptotes

### Rational [2/1] Hermite

**With f and f' at 2 points** (4 conditions):
- Numerator degree 2, denominator degree 1
- 4 DOF (3+2-1 with normalization)
- Better for functions with features like:
  - Exponential growth/decay
  - Rational functions
  - Functions approaching asymptotes

### Benchmark Results

From `special_function_convergence.py` for **exponential e^x**:

| Intervals | Cubic Hermite | [2/1] Rational | Improvement |
|-----------|---------------|----------------|-------------|
| 4         | 1.507e-03    | 6.447e-05      | **23×** |
| 8         | 8.831e-05    | 4.094e-06      | **22×** |
| 16        | 4.425e-06    | 2.569e-07      | **17×** |
| 32        | 2.104e-07    | 1.607e-08      | **13×** |

The rational approximation **consistently outperforms** polynomial for exponential functions!

## Implementation in Rust

### Direct API

```rust
use gelfgren_core::pade::TwoPointPade;

// Create [2/1] rational matching f and f' at endpoints
let left_derivs = vec![1.0, 1.0];   // [f(0), f'(0)]
let right_derivs = vec![2.718, 2.718];  // [f(1), f'(1)]

let pade = TwoPointPade::from_endpoint_derivatives(
    &left_derivs,
    &right_derivs,
    2, 1,  // [n/m]
    0.0, 1.0
)?;

let value = pade.evaluate(0.5)?;
```

### Accessing the Rational Function

```rust
// Get underlying rational
let rational = pade.rational();

// Access numerator and denominator
let numerator = rational.numerator();
let denominator = rational.denominator();

// Get Bernstein coefficients
let num_coeffs = numerator.scaled_coefficients();
let den_coeffs = denominator.scaled_coefficients();
```

## Multi-Point Extension (Future)

The Hermite module has infrastructure for **n-point** Hermite interpolation (not just 2):

```rust
use gelfgren_core::hermite::{HermiteData, LagrangeHermiteInterpolant};

// Define data at multiple points
let points = vec![0.0, 0.5, 1.0];
let values = vec![
    vec![1.0, 1.0],       // [f(0), f'(0)]
    vec![1.649, 1.649],   // [f(0.5), f'(0.5)]
    vec![2.718, 2.718],   // [f(1), f'(1)]
];

let data = HermiteData::new(points, values)?;
let interpolant = LagrangeHermiteInterpolant::new(data)?;
```

**Note**: This currently builds **polynomial** interpolants. Extending to **rational** multi-point Hermite would require implementing Traub's full n-point formulas.

## Key Takeaways

✅ **Yes**, Gelfgren supports Hermite interpolation for rational approximants
✅ **Two-point Padé** is exactly rational Hermite interpolation
✅ Can match **function values** and **arbitrary derivatives** (f', f'', f''', ...)
✅ More derivatives → higher p → better approximation
✅ Rational outperforms polynomial for many function classes
✅ **Piecewise** extension allows global approximation over multiple intervals

The constraint is always: **n + m + 1 = 2p** (must be even)

This is a powerful generalization of classical Hermite interpolation! 🎯
