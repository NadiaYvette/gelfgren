# Boundary Value Problem (BVP) Solver Status

## Current Status: Infrastructure Only ⚠️

The BVP module (`gelfgren-core/src/bvp/`) currently contains **infrastructure and design** but **no working solver implementation**.

### What Exists

#### 1. Problem Definition (`problem.rs`)
- `BoundaryValueProblem<T>` - Generic BVP specification
- `DifferentialOperator` trait - For defining L[u]
- `SecondOrderLinear` - Second-order linear ODE: u'' + p(x)u' + q(x)u = f(x)
- `RightHandSide` trait - For defining f(x)
- `FunctionRHS` - RHS from a function

#### 2. Boundary Conditions (`boundary.rs`)
- `BoundaryCondition` - Dirichlet, Neumann, Robin conditions
- Flexible specification at arbitrary points
- Support for: u(a)=α, u'(a)=α, c₀u(a)+c₁u'(a)=α

#### 3. Solver Interface (`solver.rs`)
- `BVPSolver::solve()` - **Returns NotImplemented error**
- `BVPSolver::validate()` - Validates a solution (works if solution provided)

### What's Missing: The Actual Solver!

The `solve()` method is a placeholder:

```rust
pub fn solve<T>(
    _bvp: &BoundaryValueProblem<T>,
    _mesh: Mesh<T>,
    _n: usize,
    _m: usize,
) -> Result<PiecewiseRational<T>>
{
    Err(GelfgrenError::NotImplemented(
        "Full BVP collocation solver not yet implemented".to_string(),
    ))
}
```

## Why BVP Benchmarks Don't Work

The Python BVP benchmark (`benchmarks/python/bvp_convergence.py`) can't use Gelfgren's BVP solver because:

1. **No solver exists** - `BVPSolver::solve()` returns error
2. **Python uses finite differences** - Falls back to basic polynomial finite difference method
3. **"Rational" method is fake** - Just interpolates the polynomial solution

This is why the benchmark results were identical - both methods used the same polynomial finite difference solver!

## What Would Be Needed for a Working BVP Solver

### Option 1: Collocation Method (Recommended)

Implement a collocation-based solver using piecewise rational approximation:

#### Algorithm
1. **Discretize domain**: Choose mesh τ = {x₀, x₁, ..., xₙ}
2. **Represent solution**: u(x) ≈ S_{n,m}(x) (piecewise rational)
3. **Set up linear system**:
   ```
   For interior points xᵢ:
     L[S](xᵢ) = f(xᵢ)  (enforce ODE)

   For boundary points:
     BC[S] = boundary values  (enforce BCs)
   ```
4. **Solve system**: Get derivative values at mesh points
5. **Construct solution**: Build `PiecewiseRational` from derivative values

#### Implementation Steps

**Step 1: Collocation Matrix Assembly**
```rust
fn assemble_collocation_matrix<T>(
    operator: &dyn DifferentialOperator<T>,
    mesh: &Mesh<T>,
    n: usize,
    m: usize,
) -> (Vec<T>, Vec<T>) {
    // Build matrix A and RHS b such that Ax = b
    // where x contains function/derivative values at mesh points
}
```

**Step 2: Differentiation Matrices**
- Need to compute how S'(xᵢ), S''(xᵢ) depend on coefficients
- For two-point Padé: can differentiate symbolically
- Build stencils for each collocation point

**Step 3: Boundary Condition Rows**
```rust
fn add_boundary_conditions<T>(
    matrix: &mut Vec<T>,
    rhs: &mut Vec<T>,
    bcs: &[BoundaryCondition<T>],
    mesh: &Mesh<T>,
) {
    // Add rows to enforce u(a)=α, u'(b)=β, etc.
}
```

**Step 4: Linear System Solve**
- Already have `linear_system::solve_linear_system()` ✅
- Handle potential ill-conditioning with pivoting

**Step 5: Solution Construction**
```rust
fn construct_solution<T>(
    derivative_values: &[T],
    mesh: Mesh<T>,
    n: usize,
    m: usize,
) -> Result<PiecewiseRational<T>> {
    // Build PiecewiseRational from solved derivative values
    PiecewiseRational::from_mesh(mesh, n, m)
}
```

### Option 2: Shooting Method

Convert BVP to initial value problems:
1. Guess initial slope u'(a)
2. Solve IVP forward using piecewise rational
3. Check if u(b) matches boundary condition
4. Iterate with Newton's method

**Pros:** Reuses IVP solvers
**Cons:** More complex iteration, stability issues

### Option 3: Galerkin/Finite Element

Use weak formulation with rational basis functions:
1. Multiply ODE by test function, integrate
2. Use piecewise rationals as basis
3. Get linear system for coefficients

**Pros:** Better for conservation laws
**Cons:** More complex integration

## Recommended Implementation Plan

### Phase 1: Simple Second-Order Linear BVP (3-4 days)

Focus on: u'' + p(x)u' + q(x)u = f(x) with Dirichlet BCs

1. **Implement differentiation for two-point Padé** (1 day)
   - Compute S'(xᵢ), S''(xᵢ) in terms of endpoint derivatives
   - Create differentiation matrices

2. **Assemble collocation system** (1 day)
   - Interior rows: u''(xᵢ) + p(xᵢ)u'(xᵢ) + q(xᵢ)u(xᵢ) = f(xᵢ)
   - Boundary rows: u(a) = α, u(b) = β
   - Build matrix and RHS vector

3. **Solve and construct solution** (1 day)
   - Call `solve_linear_system()`
   - Extract derivative values from solution
   - Populate mesh with values
   - Build `PiecewiseRational::from_mesh()`

4. **Test with canonical problems** (1 day)
   - u'' + u = 0, u(0)=0, u(π)=0 → solution u=0
   - u'' = -π² sin(πx), u(0)=0, u(1)=0 → solution u=sin(πx)
   - Compare against analytical solutions

### Phase 2: More General Cases (2-3 days)

5. **Support Neumann/Robin BCs**
6. **Handle variable coefficients p(x), q(x)**
7. **Add adaptive mesh refinement**
8. **Error estimation**

### Phase 3: Nonlinear Extension (1-2 weeks)

9. **Newton iteration** for nonlinear ODEs
10. **Continuation methods** for parameter studies

## Current Workaround: Special Functions Benchmark

Since BVP solver doesn't exist, we focus on **special function approximation**:
- ✅ Direct interpolation/approximation of known functions
- ✅ Tests piecewise rational vs polynomial
- ✅ Shows clear advantages (23× improvement for e^x!)

This is actually the **appropriate use case** for Gelfgren as an interpolation library.

## Alternative: Just Use Existing BVP Solvers

For actual BVP solving, users can:
1. **Use established solvers**:
   - SciPy's `solve_bvp` (Python)
   - bvp4c (MATLAB)
   - BVPSuite (Fortran)

2. **Use Gelfgren for post-processing**:
   - Solve BVP with finite differences
   - Interpolate solution with piecewise rational
   - Get smoother, more accurate representation

This is what the current Python benchmark accidentally does!

## Conclusion

### Current Reality
- ✅ BVP infrastructure (types, traits, validation)
- ❌ No actual BVP solver implementation
- ❌ BVP benchmarks use placeholder code
- ✅ Special functions benchmarks work great!

### Recommendation

**Do NOT implement BVP solver** unless specifically needed for your use case because:

1. **Complex implementation** - 1-2 weeks of work minimum
2. **Existing solutions** - Many good BVP solvers already exist
3. **Not core mission** - Gelfgren is primarily an **interpolation/approximation library**
4. **Special functions work** - Already have working benchmarks showing rational advantages

### If You Do Need BVP Solver

Follow the Phase 1 plan above for second-order linear problems. Start with:
- `u'' = f(x)` with Dirichlet BCs
- Test on Poisson equation
- Verify against analytical solutions
- Then extend to more general cases

The infrastructure is there, just needs the collocation matrix assembly and system solving code!
