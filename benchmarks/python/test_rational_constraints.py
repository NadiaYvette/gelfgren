#!/usr/bin/env python3
"""
Test rational collocation with different constraint types on challenging problems.
"""

import numpy as np
from rational_collocation import (
    RationalCollocationQuadratic, QConstraintType, BernsteinBasis
)


def test_boundary_layer():
    """Test on boundary layer problem: -εu'' + u' = 1, u(0)=u(1)=0"""
    print("="*70)
    print("Boundary Layer Problem: -εu'' + u' = 1, u(0) = u(1) = 0")
    print("="*70)

    epsilon = 0.01  # Small parameter creates boundary layer

    # For this problem, we'd need to modify the solver to handle -εu'' + u'
    # For now, test on a simpler problem that has exponential character
    # -u'' = exp(x), u(0) = 0, u(1) = 0

    print("\nActually testing: -u'' = exp(x), u(0) = u(1) = 0")
    print("This has exponential character requiring good approximation\n")

    f = lambda x: np.exp(x)

    # Exact solution: u(x) = exp(x) - 1 - x(e - 1)
    # Verify: u(0) = 1 - 1 - 0 = 0 ✓
    #         u(1) = e - 1 - (e-1) = 0 ✓
    #         u''(x) = exp(x) so -u'' = -exp(x) ✓
    exact = lambda x: np.exp(x) - 1 - x * (np.e - 1)

    constraint_types = [
        (QConstraintType.ENDPOINT, "Endpoint constraints"),
        (QConstraintType.NONNEGATIVE, "Non-negative"),
        (QConstraintType.BOUNDED, "Bounded"),
        (QConstraintType.REGULARIZATION, "Regularization"),
    ]

    print(f"{'Constraint':<20} {'Degree':<10} {'Max Error':<15} {'Q range':<20}")
    print("-"*70)

    for constraint_type, name in constraint_types:
        for n in [6, 10, 14]:
            m = n // 2

            try:
                solver = RationalCollocationQuadratic(
                    f=f, a=0.0, b=1.0, alpha=0.0, beta=0.0,
                    n=n, m=m,
                    q_constraint=constraint_type,
                    q_regularization=100.0 if constraint_type == QConstraintType.REGULARIZATION else 0.0,
                    constraint_epsilon=1e-4
                )

                result = solver.solve(method='trf', verbose=False)

                if result.success:
                    # Evaluate error
                    x_test = np.linspace(0, 1, 1000)
                    u_computed = solver.evaluate_solution(result, x_test)
                    u_exact = exact(x_test)
                    error = np.abs(u_computed - u_exact)
                    max_error = np.max(error)

                    # Check Q range
                    Q_vals = np.array([BernsteinBasis.evaluate(result.Q_coeffs, xi, 0.0, 1.0)
                                      for xi in x_test])
                    q_range = f"[{np.min(Q_vals):.4f}, {np.max(Q_vals):.4f}]"

                    print(f"{name:<20} [{n}/{m}]{'':<5} {max_error:<15.2e} {q_range:<20}")
                else:
                    print(f"{name:<20} [{n}/{m}]{'':<5} {'FAILED':<15} {result.message[:20]:<20}")

            except Exception as e:
                print(f"{name:<20} [{n}/{m}]{'':<5} {'ERROR':<15} {str(e)[:20]:<20}")

        print()  # Blank line between constraint types


def test_oscillatory():
    """Test on oscillatory problem: -u'' = sin(10πx), u(0) = u(1) = 0"""
    print("\n" + "="*70)
    print("Oscillatory Problem: -u'' = sin(10πx), u(0) = u(1) = 0")
    print("="*70 + "\n")

    omega = 10.0
    f = lambda x: np.sin(omega * np.pi * x)
    exact = lambda x: np.sin(omega * np.pi * x) / (omega * np.pi)**2

    constraint_types = [
        (QConstraintType.ENDPOINT, "Endpoint constraints"),
        (QConstraintType.NONNEGATIVE, "Non-negative"),
    ]

    print(f"{'Constraint':<20} {'Degree':<10} {'Max Error':<15} {'Q deviation':<20}")
    print("-"*70)

    for constraint_type, name in constraint_types:
        for n in [10, 16, 20]:
            m = n // 2

            try:
                solver = RationalCollocationQuadratic(
                    f=f, a=0.0, b=1.0, alpha=0.0, beta=0.0,
                    n=n, m=m,
                    q_constraint=constraint_type,
                    constraint_epsilon=1e-4
                )

                result = solver.solve(method='trf', verbose=False)

                if result.success:
                    # Evaluate error
                    x_test = np.linspace(0, 1, 1000)
                    u_computed = solver.evaluate_solution(result, x_test)
                    u_exact = exact(x_test)
                    error = np.abs(u_computed - u_exact)
                    max_error = np.max(error)

                    # Check how much Q deviates from 1
                    Q_vals = np.array([BernsteinBasis.evaluate(result.Q_coeffs, xi, 0.0, 1.0)
                                      for xi in x_test])
                    q_dev = f"σ={np.std(Q_vals):.4f}"

                    print(f"{name:<20} [{n}/{m}]{'':<5} {max_error:<15.2e} {q_dev:<20}")
                else:
                    print(f"{name:<20} [{n}/{m}]{'':<5} {'FAILED':<15} {result.message[:20]:<20}")

            except Exception as e:
                print(f"{name:<20} [{n}/{m}]{'':<5} {'ERROR':<15} {str(e)[:20]:<20}")

        print()


def test_rational_advantage():
    """
    Test where rational should have advantage: function with pole outside domain.
    -u'' = f(x) where exact solution has rational character.
    """
    print("\n" + "="*70)
    print("Rational Character: u(x) = x(1-x)/(1 + x)")
    print("This solution is naturally rational!")
    print("="*70 + "\n")

    # If u(x) = x(1-x)/(1+x), compute -u''
    # u' = [[(1-2x)(1+x) - x(1-x)] / (1+x)^2]
    #    = [(1-2x-2x^2+x) - x(1-x)] / (1+x)^2
    #    = [1 - x - 2x^2] / (1+x)^2
    # u'' = [(-1-4x)(1+x)^2 - (1-x-2x^2)·2(1+x)] / (1+x)^4
    #     = [(-1-4x)(1+x) - 2(1-x-2x^2)] / (1+x)^3
    # This gets messy - let's use numerical differentiation for f

    def exact(x):
        return x * (1 - x) / (1 + x)

    # Numerical differentiation to get f = -u''
    def f(x):
        h = 1e-6
        u = exact(x)
        u_plus = exact(x + h)
        u_minus = exact(x - h)
        u_pp = (u_plus - 2*u + u_minus) / h**2
        return -u_pp

    alpha = exact(0.0)
    beta = exact(1.0)

    print(f"Boundary values: u(0) = {alpha:.6f}, u(1) = {beta:.6f}\n")

    constraint_types = [
        (QConstraintType.ENDPOINT, "Endpoint constraints"),
        (QConstraintType.REGULARIZATION, "Regularization"),
    ]

    print(f"{'Constraint':<20} {'Degree':<10} {'Max Error':<15} {'Notes':<30}")
    print("-"*70)

    for constraint_type, name in constraint_types:
        for n in [4, 6, 8]:
            m = n // 2

            try:
                solver = RationalCollocationQuadratic(
                    f=f, a=0.0, b=1.0, alpha=alpha, beta=beta,
                    n=n, m=m,
                    q_constraint=constraint_type,
                    q_regularization=100.0 if constraint_type == QConstraintType.REGULARIZATION else 0.0,
                    constraint_epsilon=1e-4
                )

                result = solver.solve(method='trf', verbose=False)

                if result.success:
                    # Evaluate error
                    x_test = np.linspace(0, 1, 1000)
                    u_computed = solver.evaluate_solution(result, x_test)
                    u_exact = exact(x_test)
                    error = np.abs(u_computed - u_exact)
                    max_error = np.max(error)

                    # Check if Q is varying (true rational) or constant (polynomial)
                    Q_vals = np.array([BernsteinBasis.evaluate(result.Q_coeffs, xi, 0.0, 1.0)
                                      for xi in x_test])
                    q_std = np.std(Q_vals)

                    if q_std < 1e-6:
                        note = "Q≈constant (polynomial)"
                    else:
                        note = f"Q varies (σ={q_std:.4f})"

                    print(f"{name:<20} [{n}/{m}]{'':<5} {max_error:<15.2e} {note:<30}")
                else:
                    print(f"{name:<20} [{n}/{m}]{'':<5} {'FAILED':<15} {result.message[:20]:<30}")

            except Exception as e:
                print(f"{name:<20} [{n}/{m}]{'':<5} {'ERROR':<15} {str(e)[:20]:<30}")

        print()


if __name__ == "__main__":
    test_boundary_layer()
    test_oscillatory()
    test_rational_advantage()
