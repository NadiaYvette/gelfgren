# Constraint Interfaces for Optimization and Root-Finding

## Current Status: Not Implemented ❌

Currently, the Padé approximant construction **builds and solves constraints internally** without exposing them to the user. This document outlines what exists, what's missing, and how to implement constraint interfaces for optimization and root-finding.

## What Currently Exists

### 1. BVP Residual Evaluation (✅ Available)

For boundary value problems, there's a `residual()` function that evaluates how well a candidate solution satisfies the ODE:

```rust
use gelfgren_core::bvp::BoundaryValueProblem;

let bvp = /* ... define BVP ... */;
let solution = /* ... some PiecewiseRational ... */;

// Evaluate at a point
let x = 0.5;
let u = solution.evaluate(x)?;
let u_prime = /* compute derivative */;
let u_double_prime = /* compute second derivative */;

// Residual: L[u](x) - f(x)
let residual = bvp.residual(x, u, &[u_prime, u_double_prime]);
```

This is useful for **validating** a solution but doesn't help **construct** one.

### 2. Internal Linear System (❌ Not Exposed)

The two-point Padé construction internally builds a linear system:

```rust
// Inside construct_rational_pade():
let matrix_a: Vec<T> = /* ... build constraint matrix ... */;
let vector_b: Vec<T> = /* ... build RHS ... */;

// Immediately solve
let solution = solve_linear_system(&matrix_a, &vector_b, num_unknowns)?;
```

But these constraints are **not exposed** to the user!

## What's Missing: Constraint Function Interface

To enable optimization and root-finding, we need interfaces like:

### Desired API 1: Residual Function

```rust
pub struct HermiteConstraints<T> {
    // Configuration
    left_derivatives: Vec<T>,
    right_derivatives: Vec<T>,
    n: usize,  // numerator degree
    m: usize,  // denominator degree
    x0: T,
    x1: T,
}

impl<T: Float> HermiteConstraints<T> {
    /// Evaluate constraint residuals given coefficient vector.
    ///
    /// # Arguments
    /// * `coeffs` - [a₀, ..., aₙ, b₁, ..., bₘ] (b₀ = 1 - Σbᵢ from normalization)
    ///
    /// # Returns
    /// Vector of residuals r where ideally r = 0
    pub fn residuals(&self, coeffs: &[T]) -> Vec<T> {
        // For each matching condition:
        // 1. Compute R^(k)(xᵢ) from coefficients
        // 2. Compare to target f^(k)(xᵢ)
        // 3. Return difference
    }

    /// Objective function: sum of squared residuals.
    pub fn objective(&self, coeffs: &[T]) -> T {
        let residuals = self.residuals(coeffs);
        residuals.iter().map(|r| r * r).sum()
    }

    /// Gradient of objective (for optimization).
    pub fn gradient(&self, coeffs: &[T]) -> Vec<T> {
        // ∂/∂cᵢ Σⱼ rⱼ²
    }

    /// Constraint Jacobian (for root-finding).
    pub fn jacobian(&self, coeffs: &[T]) -> Vec<Vec<T>> {
        // ∂rᵢ/∂cⱼ for all i,j
    }
}
```

### Desired API 2: Linear System Extraction

```rust
impl<T: Float> HermiteConstraints<T> {
    /// Build the linear system Ax = b without solving it.
    pub fn build_linear_system(&self) -> (Vec<T>, Vec<T>, usize) {
        // Returns (matrix_a, vector_b, num_unknowns)
        // User can solve with their own method or inspect structure
    }

    /// Check if system is well-conditioned.
    pub fn condition_number(&self) -> T {
        let (a, _, n) = self.build_linear_system();
        // Compute κ(A) = ||A|| ||A⁻¹||
    }
}
```

### Desired API 3: Quadratic Form (for Optimization)

For least-squares problems, express as:
```
minimize: ||Ax - b||²
```

This is already quadratic in x!

```rust
impl<T: Float> HermiteConstraints<T> {
    /// Get quadratic form matrices: minimize x^T H x + g^T x
    /// where H = A^T A, g = -2 A^T b
    pub fn quadratic_form(&self) -> (Vec<T>, Vec<T>) {
        let (a, b, n) = self.build_linear_system();

        // H = A^T A (Hessian)
        let h = matrix_multiply_transpose(a, a, n, n);

        // g = -2 A^T b (gradient)
        let g = matrix_vector_multiply_transpose(a, b, n);
        let g: Vec<_> = g.iter().map(|gi| *gi * T::from(-2.0).unwrap()).collect();

        (h, g)
    }
}
```

## Use Cases

### 1. Root-Finding for Nonlinear Hermite Interpolation

When the matching conditions are nonlinear (e.g., interpolating at points determined by the solution itself):

```rust
use gelfgren::optimize::newton_raphson;

let constraints = HermiteConstraints::new(/* ... */);

// Initial guess for coefficients
let x0 = vec![/* initial guess */];

// Solve F(x) = 0 using Newton's method
let solution = newton_raphson(
    |x| constraints.residuals(x),
    |x| constraints.jacobian(x),
    x0,
    max_iter = 100,
    tol = 1e-10
)?;

// Build rational from solution
let rational = constraints.build_from_coefficients(&solution)?;
```

### 2. Constrained Optimization

Minimize an objective while satisfying Hermite constraints:

```rust
use gelfgren::optimize::constrained_minimize;

// Objective: minimize ∫[0,1] |R(x) - target(x)|² dx
let objective = |coeffs: &[f64]| {
    let rational = build_rational_from_coeffs(coeffs);
    integrate(|x| {
        let r = rational.evaluate(x);
        let t = target(x);
        (r - t) * (r - t)
    }, 0.0, 1.0)
};

// Constraints: Hermite matching conditions
let constraints = HermiteConstraints::new(/* ... */);

let solution = constrained_minimize(
    objective,
    |x| constraints.residuals(x),  // equality constraints
    x0,
)?;
```

### 3. Overdetermined Systems (Least Squares)

When you have more constraints than degrees of freedom:

```rust
// Example: 6 constraints but only [2/1] rational (4 DOF)
let left_derivs = vec![f(0), f_prime(0), f_double_prime(0)];  // 3 values
let right_derivs = vec![f(1), f_prime(1), f_double_prime(1)]; // 3 values

let constraints = HermiteConstraints::new(
    left_derivs, right_derivs,
    2, 1,  // [2/1] has only 4 DOF, but we have 6 constraints!
    0.0, 1.0
);

// Cannot satisfy exactly, use least squares:
// minimize ||Ax - b||²
let (h, g) = constraints.quadratic_form();
let solution = quadratic_minimize(h, g)?;

// Get "best fit" rational that approximately satisfies constraints
let rational = constraints.build_from_coefficients(&solution)?;
```

### 4. Underdetermined Systems (Multiple Solutions)

When you have fewer constraints than DOF:

```rust
// Example: Only match function values (not derivatives)
let left_derivs = vec![f(0)];   // 1 constraint
let right_derivs = vec![f(1)];  // 1 constraint

// [2/1] has 4 DOF, but only 2 constraints!

let constraints = HermiteConstraints::new(
    left_derivs, right_derivs,
    2, 1,
    0.0, 1.0
);

// Add regularization to pick smooth solution:
// minimize ||R''||² subject to constraints
let objective = |coeffs: &[f64]| {
    let rational = build_rational_from_coeffs(coeffs);
    integrate(|x| {
        let r_double_prime = rational.derivative_nth(x, 2);
        r_double_prime * r_double_prime
    }, 0.0, 1.0)
};

let solution = constrained_minimize(
    objective,
    |x| constraints.residuals(x),
    x0,
)?;
```

### 5. BVP as Optimization Problem

Transform BVP to optimization:

```rust
// BVP: u'' + u = 0, u(0)=0, u(π)=0
let bvp = BoundaryValueProblem::new(/* ... */);
let mesh = Mesh::uniform(0.0, PI, 10);

// Objective: minimize PDE residual + BC violation
let objective = |coeffs: &[f64]| {
    let solution = build_piecewise_rational_from_coeffs(coeffs, &mesh, 2, 1);

    // PDE residual
    let mut pde_error = 0.0;
    for x in mesh.interior_points() {
        let residual = bvp.residual(x, solution.evaluate(x), /* ... */);
        pde_error += residual * residual;
    }

    // Boundary condition violation
    let bc_error = (solution.evaluate(0.0) - 0.0).powi(2)
                 + (solution.evaluate(PI) - 0.0).powi(2);

    pde_error + 1000.0 * bc_error  // Weight BC heavily
};

// No explicit constraints, pure optimization
let solution = unconstrained_minimize(objective, x0)?;
```

## Implementation Plan

### Phase 1: Extract Linear System (1 day)

Refactor `construct_rational_pade` to expose the linear system:

```rust
// Current (internal):
fn construct_rational_pade(/* ... */) -> Result<(Poly, Poly)> {
    // Build system
    let (a, b, n) = build_linear_system(/* ... */);
    // Solve immediately
    let solution = solve_linear_system(&a, &b, n)?;
    // Return polynomials
}

// New (exposed):
pub fn build_hermite_linear_system<T>(
    left_derivatives: &[T],
    right_derivatives: &[T],
    n: usize,
    m: usize,
    x0: T,
    x1: T,
) -> Result<(Vec<T>, Vec<T>, usize)> {
    // Just build and return, don't solve
}

// Use in construct_rational_pade:
fn construct_rational_pade(/* ... */) -> Result<(Poly, Poly)> {
    let (a, b, num) = build_hermite_linear_system(/* ... */)?;
    let solution = solve_linear_system(&a, &b, num)?;
    // Return polynomials
}
```

### Phase 2: Residual Function (2 days)

Implement `HermiteConstraints` struct with residual evaluation:

```rust
pub struct HermiteConstraints<T> { /* ... */ }

impl<T> HermiteConstraints<T> {
    pub fn new(/* ... */) -> Result<Self>;

    pub fn residuals(&self, coeffs: &[T]) -> Vec<T> {
        // 1. Unpack coefficients into numerator/denominator
        // 2. Build Bernstein polynomials
        // 3. Build rational R(x) = P(x)/Q(x)
        // 4. For each constraint k at each endpoint i:
        //    - Compute R^(k)(xᵢ)
        //    - Compare to target f^(k)(xᵢ)
        //    - Store residual
        // 5. Return residual vector
    }
}
```

### Phase 3: Jacobian (3 days)

Implement analytical Jacobian for Newton methods:

```rust
impl<T> HermiteConstraints<T> {
    pub fn jacobian(&self, coeffs: &[T]) -> Vec<Vec<T>> {
        // For each residual rᵢ:
        //   For each coefficient cⱼ:
        //     Compute ∂rᵢ/∂cⱼ analytically
        //
        // This requires derivatives of R^(k)(x) with respect to coefficients
    }
}
```

### Phase 4: Quadratic Form (1 day)

For least-squares optimization:

```rust
impl<T> HermiteConstraints<T> {
    pub fn quadratic_form(&self) -> (Vec<T>, Vec<T>) {
        let (a, b, n) = self.build_linear_system();
        // H = A^T A
        // g = -2 A^T b
    }
}
```

### Phase 5: Integration with Optimizers (2 days)

Create adapter functions for common optimization libraries:

```rust
// For scipy.optimize (Python)
#[pyclass]
struct PyHermiteConstraints {
    inner: HermiteConstraints<f64>,
}

#[pymethods]
impl PyHermiteConstraints {
    fn residuals(&self, coeffs: PyReadonlyArray1<f64>) -> Py<PyArray1<f64>> {
        // Call Rust, return to Python
    }

    fn jacobian(&self, coeffs: PyReadonlyArray1<f64>) -> Py<PyArray2<f64>> {
        // Call Rust, return to Python
    }
}

// Python usage:
// from scipy.optimize import least_squares
// constraints = gf.HermiteConstraints(left, right, n, m, a, b)
// result = least_squares(constraints.residuals, x0, jac=constraints.jacobian)
```

## Benefits

### 1. **Flexibility**
- Users can choose their own solvers (Newton, trust-region, BFGS, etc.)
- Handle special cases (sparse, structured, etc.)

### 2. **Nonlinear Extensions**
- Currently limited to linear matching conditions
- With constraint interface, can handle nonlinear problems

### 3. **Constrained Problems**
- Add additional constraints beyond Hermite matching
- E.g., R(x) ≥ 0, ∫R(x)dx = 1, etc.

### 4. **Analysis**
- Inspect condition number before solving
- Understand constraint structure
- Debug ill-conditioned problems

### 5. **BVP Optimization**
- Transform BVPs to optimization problems
- Use modern optimization methods
- Handle nonlinear BVPs

## Recommended Priority

1. **Phase 1** (HIGH): Extract linear system - enables advanced users immediately
2. **Phase 2** (HIGH): Residual function - core functionality for root-finding
3. **Phase 4** (MEDIUM): Quadratic form - useful for least-squares
4. **Phase 3** (LOW): Jacobian - optimization, but finite differences work
5. **Phase 5** (LOW): Python integration - nice to have

## Current Workaround

For now, users who need constraint interfaces can:

1. **Use finite differences** for gradients/Jacobians
2. **Build their own linear system** by studying the code
3. **Use existing solver** and post-process if needed

But this is clearly **suboptimal** - a proper constraint interface would be very valuable!

## Conclusion

**Status**: ❌ Not currently implemented
**Priority**: 🔥 HIGH for advanced use cases
**Effort**: ~1-2 weeks for full implementation
**Value**: Enables optimization, root-finding, constrained problems, BVP solving via optimization

This would be a **significant enhancement** making Gelfgren useful for a much broader class of problems! 🚀
