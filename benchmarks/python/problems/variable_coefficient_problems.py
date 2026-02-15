#!/usr/bin/env python3
"""
Variable Coefficient BVP Test Problems

Test problems for: -(a(x)*u')' = f(x) on [0, 1]
                   u(0) = α, u(1) = β

These problems test method performance on variable diffusion coefficients.
CRITICAL: Discontinuous a(x) should cause catastrophic failure for rationals.
"""

import numpy as np
from problem_types import create_variable_coefficient_problem


def variable_coefficient_polynomial():
    """
    VC1: Smooth polynomial coefficient

    a(x) = 1 + x²
    Exact: u(x) = sin(πx)

    Smooth varying coefficient. All methods should work well.
    """
    def a_func(x):
        return 1 + x**2

    def exact(x):
        return np.sin(np.pi * x)

    def forcing(x):
        # u = sin(πx)
        # u' = π*cos(πx)
        # a(x)*u' = (1 + x²)*π*cos(πx)
        # (a(x)*u')' = π*[-sin(πx)] * (1 + x²) + π*cos(πx)*2x
        #            = -π*(1 + x²)*sin(πx) + 2πx*cos(πx)
        # -(a(x)*u')' = π*(1 + x²)*sin(πx) - 2πx*cos(πx)
        return np.pi * (1 + x**2) * np.sin(np.pi * x) - 2 * np.pi * x * np.cos(np.pi * x)

    return create_variable_coefficient_problem(
        name='VarCoeff_Poly',
        domain=(0.0, 1.0),
        a_func=a_func,
        forcing=forcing,
        exact=exact,
        alpha=0.0,
        beta=0.0,
        description='Smooth polynomial coefficient a(x)=1+x², all methods should work well',
        expected_difficulty='smooth'
    )


def variable_coefficient_exponential():
    """
    VC2: Smooth exponential coefficient

    a(x) = e^x
    Exact: u(x) = x(1-x)

    Smooth exponentially varying coefficient.
    Rationals should handle natural exponential representation.
    """
    def a_func(x):
        return np.exp(x)

    def exact(x):
        return x * (1 - x)

    def forcing(x):
        # u = x - x²
        # u' = 1 - 2x
        # a(x)*u' = e^x * (1 - 2x)
        # (a(x)*u')' = e^x*(1 - 2x) + e^x*(-2) = e^x*(1 - 2x - 2) = e^x*(-1 - 2x)
        # -(a(x)*u')' = e^x*(1 + 2x)
        return np.exp(x) * (1 + 2 * x)

    return create_variable_coefficient_problem(
        name='VarCoeff_Exp',
        domain=(0.0, 1.0),
        a_func=a_func,
        forcing=forcing,
        exact=exact,
        alpha=0.0,
        beta=0.0,
        description='Smooth exponential coefficient a(x)=e^x, rationals should work well',
        expected_difficulty='smooth'
    )


def variable_coefficient_discontinuous():
    """
    VC3: Discontinuous coefficient (jump at x=0.5)

    a(x) = 1 if x < 0.5, else 10
    Exact solution manufactured to have continuous u and u' but discontinuous u''.

    **CRITICAL: Rationals expected to CATASTROPHICALLY FAIL**

    Rational approximants cannot represent discontinuities without Gibbs phenomena.
    This is a key failure mode.
    """
    def a_func(x):
        if hasattr(x, '__len__'):
            result = np.ones_like(x)
            result[x >= 0.5] = 10.0
            return result
        else:
            return 1.0 if x < 0.5 else 10.0

    # Manufacture solution with continuous u and u', but a(x)*u' has jump
    # Let u be smooth across interface, but derivatives scaled by a(x)
    def exact(x):
        # Piecewise smooth solution
        # Left (x < 0.5): u = x²
        # Right (x >= 0.5): u = -0.25 + x (continuous at 0.5)
        # Actually this won't have continuous derivative...

        # Better: Make u and u' continuous
        # Left: u = x²
        # At x=0.5: u(0.5) = 0.25, u'(0.5) = 1
        # Right: u = ax² + bx + c with u(0.5)=0.25, u'(0.5)=1, u(1)=0
        # a + b/2 + c = 0.25
        # 2a*0.5 + b = 1 => a + b = 1
        # a + b + c = 0
        # From a+b=1: b = 1-a
        # From a+b+c=0: c = -1
        # From a+b/2+c=0.25: a + (1-a)/2 - 1 = 0.25 => a + 0.5 - a/2 - 1 = 0.25
        #                                           => a/2 = 0.75 => a = 1.5
        # So a=1.5, b=-0.5, c=-1
        # Right: u = 1.5x² - 0.5x - 1

        if hasattr(x, '__len__'):
            result = np.zeros_like(x)
            left = x < 0.5
            right = x >= 0.5
            result[left] = x[left]**2
            result[right] = 1.5 * x[right]**2 - 0.5 * x[right] - 1.0
            return result
        else:
            if x < 0.5:
                return x**2
            else:
                return 1.5 * x**2 - 0.5 * x - 1.0

    def forcing(x):
        # -(a(x)*u')' = f
        # Left (x < 0.5): u = x², u' = 2x, a*u' = 1*2x = 2x
        #                 (a*u')' = 2, -(a*u')' = -2
        # Right (x >= 0.5): u = 1.5x² - 0.5x - 1, u' = 3x - 0.5, a*u' = 10*(3x - 0.5)
        #                   (a*u')' = 10*3 = 30, -(a*u')' = -30

        if hasattr(x, '__len__'):
            result = np.zeros_like(x)
            left = x < 0.5
            right = x >= 0.5
            result[left] = -2.0
            result[right] = -30.0
            return result
        else:
            return -2.0 if x < 0.5 else -30.0

    return create_variable_coefficient_problem(
        name='VarCoeff_Discontinuous',
        domain=(0.0, 1.0),
        a_func=a_func,
        forcing=forcing,
        exact=exact,
        alpha=0.0,
        beta=0.0,
        description='Discontinuous coefficient (jump at x=0.5), **RATIONALS FAIL CATASTROPHICALLY**',
        expected_difficulty='discontinuous'
    )


def variable_coefficient_oscillatory():
    """
    VC4: Rapidly oscillating coefficient

    a(x) = 2 + sin(20πx)
    Exact: u(x) = x(1-x)

    Rapidly varying smooth coefficient (10 oscillations).
    Tests resolution requirements; all methods need fine grids.
    """
    def a_func(x):
        return 2 + np.sin(20 * np.pi * x)

    def exact(x):
        return x * (1 - x)

    def forcing(x):
        # u = x - x²
        # u' = 1 - 2x
        # a(x) = 2 + sin(20πx)
        # a'(x) = 20π*cos(20πx)
        # a(x)*u' = [2 + sin(20πx)] * (1 - 2x)
        # (a(x)*u')' = a'(x)*u' + a(x)*u''
        #            = 20π*cos(20πx)*(1 - 2x) + [2 + sin(20πx)]*(-2)
        # -(a(x)*u')' = -20π*cos(20πx)*(1 - 2x) + 2*[2 + sin(20πx)]
        return -20 * np.pi * np.cos(20 * np.pi * x) * (1 - 2 * x) + 2 * (2 + np.sin(20 * np.pi * x))

    return create_variable_coefficient_problem(
        name='VarCoeff_Oscillatory',
        domain=(0.0, 1.0),
        a_func=a_func,
        forcing=forcing,
        exact=exact,
        alpha=0.0,
        beta=0.0,
        description='Rapidly oscillating coefficient a(x)=2+sin(20πx), high resolution required',
        expected_difficulty='oscillatory'
    )


def variable_coefficient_near_singular():
    """
    VC5: Near-singular coefficient

    a(x) = √(x + 0.01)
    Exact: u(x) = x²

    Coefficient nearly singular at x=0 (derivative blows up).
    Tests method robustness near singularities.
    """
    def a_func(x):
        return np.sqrt(x + 0.01)

    def exact(x):
        return x**2

    def forcing(x):
        # u = x²
        # u' = 2x
        # a(x) = √(x + 0.01)
        # a'(x) = 1/(2√(x + 0.01))
        # a(x)*u' = √(x + 0.01) * 2x = 2x√(x + 0.01)
        # (a(x)*u')' = a'(x)*u' + a(x)*u''
        #            = [1/(2√(x + 0.01))] * 2x + √(x + 0.01) * 2
        #            = x/√(x + 0.01) + 2√(x + 0.01)
        # -(a(x)*u')' = -x/√(x + 0.01) - 2√(x + 0.01)
        sqrt_term = np.sqrt(x + 0.01)
        return -x / sqrt_term - 2 * sqrt_term

    return create_variable_coefficient_problem(
        name='VarCoeff_NearSingular',
        domain=(0.0, 1.0),
        a_func=a_func,
        forcing=forcing,
        exact=exact,
        alpha=0.0,
        beta=1.0,
        description='Near-singular coefficient a(x)=√(x+0.01), tests robustness',
        expected_difficulty='singular'
    )


def get_all_variable_coefficient_problems():
    """Return all variable coefficient test problems."""
    return [
        variable_coefficient_polynomial(),
        variable_coefficient_exponential(),
        variable_coefficient_discontinuous(),  # CRITICAL failure test
        variable_coefficient_oscillatory(),
        variable_coefficient_near_singular(),
    ]


if __name__ == "__main__":
    # Test problem definitions
    problems = get_all_variable_coefficient_problems()

    print("Variable Coefficient Test Problems")
    print("=" * 70)

    for prob in problems:
        print(f"\n{prob.name}:")
        print(f"  Domain: [{prob.domain[0]}, {prob.domain[1]}]")
        print(f"  BCs: u({prob.domain[0]}) = {prob.boundary_conditions['alpha']:.6f}, "
              f"u({prob.domain[1]}) = {prob.boundary_conditions['beta']:.6f}")
        print(f"  Description: {prob.description}")
        print(f"  Difficulty: {prob.expected_difficulty}")

        # Test that exact solution satisfies BCs
        u_a = prob.exact(prob.domain[0])
        u_b = prob.exact(prob.domain[1])
        print(f"  Verification: u(a) = {u_a:.6f} (should be {prob.boundary_conditions['alpha']:.6f})")
        print(f"                u(b) = {u_b:.6f} (should be {prob.boundary_conditions['beta']:.6f})")

        # Sample a(x) at a few points
        x_sample = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
        a_sample = prob.coeffs['a_func'](x_sample)
        print(f"  a(x) samples: {a_sample}")
