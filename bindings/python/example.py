#!/usr/bin/env python3
"""
Example usage of Gelfgren Python bindings.

Demonstrates polynomial operations, rational functions, and Padé approximants.
"""

import numpy as np
import gelfgren

def demonstrate_polynomials():
    """Demonstrate Bernstein polynomial operations."""
    print("Bernstein Polynomial Operations")
    print("=" * 40)
    print()

    # Create polynomial P(x) = 1 + 2x + 3x²
    p = gelfgren.BernsteinPolynomial([1.0, 2.0, 3.0], 0.0, 1.0)
    print(f"Created: {p}")
    print(f"Degree: {p.degree()}")
    print(f"Interval: {p.interval()}")
    print()

    # Evaluate at multiple points
    x = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    y = p.evaluate(x)

    print("Evaluations:")
    for xi, yi in zip(x, y):
        print(f"  P({xi:.2f}) = {yi:.6f}")
    print()

    # Compute derivative
    p_prime = p.derivative()
    y_prime = p_prime.evaluate(x)

    print("Derivative P'(x):")
    for xi, yi in zip(x, y_prime):
        print(f"  P'({xi:.2f}) = {yi:.6f}")
    print()

    # Polynomial arithmetic
    q = gelfgren.BernsteinPolynomial([1.0, 0.0], 0.0, 1.0)  # q(x) = 1
    sum_poly = p + q
    diff_poly = p - q
    prod_poly = p * q

    x_test = 0.5
    print(f"Polynomial arithmetic at x = {x_test}:")
    print(f"  P(x) + Q(x) = {sum_poly.eval_scalar(x_test):.6f}")
    print(f"  P(x) - Q(x) = {diff_poly.eval_scalar(x_test):.6f}")
    print(f"  P(x) * Q(x) = {prod_poly.eval_scalar(x_test):.6f}")
    print()


def demonstrate_rational_functions():
    """Demonstrate rational function operations."""
    print("\nRational Functions")
    print("=" * 40)
    print()

    # Create R(x) = (1 + x) / (1 + 2x)
    num = gelfgren.BernsteinPolynomial([1.0, 1.0], 0.0, 1.0)
    den = gelfgren.BernsteinPolynomial([1.0, 2.0], 0.0, 1.0)
    r = gelfgren.RationalFunction(num, den)

    print("Created rational function R(x) = (1 + x) / (1 + 2x)")
    print()

    x = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    y = r.evaluate(x)

    print("Evaluations:")
    for xi, yi in zip(x, y):
        print(f"  R({xi:.1f}) = {yi:.6f}")
    print()

    # Derivative
    r_prime = r.derivative()
    y_prime = r_prime.evaluate(x)

    print("Derivative R'(x):")
    for xi, yi in zip(x, y_prime):
        print(f"  R'({xi:.1f}) = {yi:.6f}")
    print()


def demonstrate_pade_approximants():
    """Demonstrate Padé approximants."""
    print("\nPadé Approximants")
    print("=" * 40)
    print()

    # Approximate exp(x) using Taylor series
    # exp(x) = 1 + x + x²/2! + x³/3! + x⁴/4! + ...
    coeffs = [1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0]

    # Create [2/2] Padé approximant
    pade = gelfgren.PadeApproximant(coeffs, 2, 2, 0.0, -1.0, 1.0)

    print("Padé [2/2] approximant for exp(x) on [-1, 1]")
    print()

    x = np.linspace(-1.0, 1.0, 9)
    exact = np.exp(x)
    approx = pade.evaluate(x)
    error = np.abs(exact - approx)

    print(f"{'x':>8} {'exp(x)':>12} {'Padé(x)':>12} {'Error':>12}")
    print("-" * 48)
    for xi, ei, ai, erri in zip(x, exact, approx, error):
        print(f"{xi:8.2f} {ei:12.6f} {ai:12.6f} {erri:12.2e}")
    print()


def demonstrate_meshes():
    """Demonstrate mesh generation."""
    print("\nMesh Generation")
    print("=" * 40)
    print()

    # Create uniform mesh
    uniform = gelfgren.Mesh.uniform(0.0, 1.0, 4)
    print(f"Uniform mesh: {uniform}")

    # Create Chebyshev mesh
    chebyshev = gelfgren.Mesh.chebyshev(0.0, 1.0, 4)
    print(f"Chebyshev mesh: {chebyshev}")
    print("  (Chebyshev meshes cluster points near boundaries)")
    print()


def demonstrate_exception_handling():
    """Demonstrate exception handling."""
    print("\nException Handling")
    print("=" * 40)
    print()

    try:
        # Create rational function with pole
        num = gelfgren.BernsteinPolynomial([1.0, 1.0], 0.0, 1.0)
        den = gelfgren.BernsteinPolynomial([1.0, -1.0], 0.0, 1.0)
        r = gelfgren.RationalFunction(num, den)

        print("Evaluating rational function with pole...")

        # This should work
        result = r.eval_scalar(0.0)
        print(f"  R(0.0) = {result:.6f} (OK)")

        # This should raise an exception (pole near x = 0.5)
        print("  R(0.5) = ", end="")
        result = r.eval_scalar(0.5)
        print(f"{result:.6f}")

    except RuntimeError as e:
        print(f"Exception caught: {e}")
        print("  (This is expected behavior for poles)")
    print()


def main():
    """Run all examples."""
    print("Gelfgren Python Bindings Example")
    print("=" * 40)
    print()

    demonstrate_polynomials()
    demonstrate_rational_functions()
    demonstrate_pade_approximants()
    demonstrate_meshes()
    demonstrate_exception_handling()

    print("All examples completed successfully!")


if __name__ == "__main__":
    main()
