#!/usr/bin/env python3
"""
Demonstration of HermiteConstraints interface for optimization and root-finding.

This example shows how to use the constraint interface to:
1. Solve linear systems directly
2. Use root-finding (scipy.optimize.least_squares)
3. Use as an optimization objective
4. Handle overdetermined systems
"""

import numpy as np

# Note: Requires building with:
# PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop
try:
    import gelfgren as gf
except ImportError:
    print("Error: gelfgren Python bindings not installed")
    print("Build with: cd bindings/python && PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin develop")
    exit(1)


def example_1_linear_system():
    """Approach 1: Solve linear system directly."""
    print("\n=== Example 1: Direct Linear System Solution ===")

    # Create constraints for e^x on [0,1] with [2/1] rational
    left = np.array([1.0, 1.0])           # [f(0), f'(0)]
    right = np.array([np.e, np.e])        # [f(1), f'(1)]

    constraints = gf.HermiteConstraints(left, right, 2, 1, 0.0, 1.0)

    print(f"Constraints: [{constraints.numerator_degree()}/{constraints.denominator_degree()}] rational")
    print(f"System: {constraints.num_constraints()} equations, {constraints.num_unknowns()} unknowns")

    # Get linear system Ax = b
    A_flat, b, n = constraints.linear_system()

    # Reshape A to 2D matrix
    A = A_flat.reshape((constraints.num_constraints(), n))

    print(f"Matrix A shape: {A.shape}")
    print(f"Vector b shape: {b.shape}")

    # Solve directly
    solution = np.linalg.solve(A, b)
    print(f"Solution coefficients: {solution}")

    # Build rational from solution
    rational = constraints.build_rational(solution)

    # Evaluate at test points
    x_test = np.array([0.0, 0.5, 1.0])
    y_test = np.array([rational.eval_scalar(x) for x in x_test])
    y_exact = np.exp(x_test)

    print(f"\nEvaluation:")
    print(f"x:      {x_test}")
    print(f"R(x):   {y_test}")
    print(f"e^x:    {y_exact}")
    print(f"Error:  {np.abs(y_test - y_exact)}")


def example_2_root_finding():
    """Approach 2: Use scipy root-finding."""
    print("\n=== Example 2: Root-Finding with least_squares ===")

    try:
        from scipy.optimize import least_squares
    except ImportError:
        print("scipy not installed, skipping this example")
        return

    left = np.array([1.0, 1.0])
    right = np.array([np.e, np.e])

    constraints = gf.HermiteConstraints(left, right, 2, 1, 0.0, 1.0)

    # Initial guess
    x0 = np.ones(4)

    # Solve using least_squares (Levenberg-Marquardt)
    result = least_squares(
        constraints.residuals,
        x0,
        jac=constraints.jacobian,  # Use analytical Jacobian
        verbose=1
    )

    print(f"\nOptimization result:")
    print(f"Success: {result.success}")
    print(f"Iterations: {result.nfev}")
    print(f"Final residual: {result.cost:.2e}")
    print(f"Solution: {result.x}")

    # Build rational
    rational = constraints.build_rational(result.x)
    x_test = np.linspace(0, 1, 5)
    y_test = np.array([rational.eval_scalar(x) for x in x_test])
    y_exact = np.exp(x_test)

    print(f"\nMax error: {np.max(np.abs(y_test - y_exact)):.2e}")


def example_3_optimization():
    """Approach 3: Use as optimization objective."""
    print("\n=== Example 3: Optimization with minimize ===")

    try:
        from scipy.optimize import minimize
    except ImportError:
        print("scipy not installed, skipping this example")
        return

    left = np.array([1.0, 1.0])
    right = np.array([np.e, np.e])

    constraints = gf.HermiteConstraints(left, right, 2, 1, 0.0, 1.0)

    # Initial guess
    x0 = np.ones(4)

    # Minimize sum of squared residuals
    result = minimize(
        constraints.objective,
        x0,
        method='BFGS',
        options={'disp': True}
    )

    print(f"\nOptimization result:")
    print(f"Success: {result.success}")
    print(f"Iterations: {result.nit}")
    print(f"Final objective: {result.fun:.2e}")
    print(f"Solution: {result.x}")


def example_4_overdetermined():
    """Approach 4: Overdetermined system (least squares)."""
    print("\n=== Example 4: Overdetermined System ===")

    # More constraints than degrees of freedom
    left = np.array([1.0, 1.0, 1.0])      # 3 derivatives at left
    right = np.array([np.e, np.e, np.e])  # 3 derivatives at right
    # Total: 6 constraints but [2/1] has only 4 DOF!

    # This should fail the validation since n+m+1 ≠ 2p
    try:
        constraints = gf.HermiteConstraints(left, right, 2, 1, 0.0, 1.0)
        print("ERROR: Should have raised an exception!")
    except ValueError as e:
        print(f"Expected error (mismatch): {e}")

    # To handle overdetermined systems, we need to adjust n and m
    # For 6 constraints (p=3), we need n+m+1=6
    # For example, [3/2] rational: 3+2+1=6
    print("\nUsing [3/2] rational to match 6 constraints:")
    constraints = gf.HermiteConstraints(left, right, 3, 2, 0.0, 1.0)

    print(f"System: {constraints.num_constraints()} equations, {constraints.num_unknowns()} unknowns")

    # Solve
    A_flat, b, n = constraints.linear_system()
    A = A_flat.reshape((constraints.num_constraints(), n))
    solution = np.linalg.solve(A, b)

    print(f"Solution: {solution}")

    # Evaluate
    rational = constraints.build_rational(solution)
    x_test = np.array([0.0, 0.5, 1.0])
    y_test = np.array([rational.eval_scalar(x) for x in x_test])
    y_exact = np.exp(x_test)

    print(f"\nEvaluation:")
    print(f"x:      {x_test}")
    print(f"R(x):   {y_test}")
    print(f"e^x:    {y_exact}")
    print(f"Error:  {np.abs(y_test - y_exact)}")


def example_5_inspection():
    """Approach 5: Inspect constraint structure."""
    print("\n=== Example 5: Constraint Inspection ===")

    left = np.array([1.0, 1.0])
    right = np.array([np.e, np.e])

    constraints = gf.HermiteConstraints(left, right, 2, 1, 0.0, 1.0)

    # Get quadratic form
    H, g = constraints.quadratic_form()

    print(f"Hessian H shape: {H.shape}")
    print(f"Gradient g shape: {g.shape}")

    # Check if Hessian is symmetric
    is_symmetric = np.allclose(H, H.T)
    print(f"Hessian is symmetric: {is_symmetric}")

    # Compute condition number
    cond = np.linalg.cond(H)
    print(f"Condition number: {cond:.2e}")

    # Get linear system
    A_flat, b, n = constraints.linear_system()
    A = A_flat.reshape((constraints.num_constraints(), n))

    # Verify H = A^T A
    H_computed = A.T @ A
    is_correct = np.allclose(H, H_computed)
    print(f"H = A^T A verified: {is_correct}")

    # Verify g = -2 A^T b
    g_computed = -2 * A.T @ b
    is_correct_g = np.allclose(g, g_computed)
    print(f"g = -2 A^T b verified: {is_correct_g}")


if __name__ == "__main__":
    print("=" * 70)
    print("Gelfgren HermiteConstraints Interface Demonstration")
    print("=" * 70)

    example_1_linear_system()
    example_2_root_finding()
    example_3_optimization()
    example_4_overdetermined()
    example_5_inspection()

    print("\n" + "=" * 70)
    print("All examples completed!")
    print("=" * 70)
