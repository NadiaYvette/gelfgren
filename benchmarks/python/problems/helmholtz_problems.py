#!/usr/bin/env python3
"""
Helmholtz BVP Test Problems

Test problems for: -u'' + k²*u = f(x) on [0, 1]
                   u(0) = α, u(1) = β

These problems test rational collocation on problems with different
frequency content and compare against spectral/high-order methods.
"""

import numpy as np
from problem_types import create_helmholtz_problem


def helmholtz_sin_k1():
    """
    H1: Low frequency sine (-u'' + u = f)

    Exact: u(x) = sin(πx)
    k² = 1 (k = 1)

    Smooth baseline problem. All methods should converge exponentially.
    """
    k = 1.0
    k_squared = k**2

    def exact(x):
        return np.sin(np.pi * x)

    def forcing(x):
        # -u'' + k²*u = f
        # -(-π²*sin(πx)) + k²*sin(πx) = f
        return (np.pi**2 + k_squared) * np.sin(np.pi * x)

    return create_helmholtz_problem(
        name='Helmholtz_Sin_k1',
        domain=(0.0, 1.0),
        k_squared=k_squared,
        forcing=forcing,
        exact=exact,
        alpha=0.0,
        beta=0.0,
        description=f'Low frequency sine, k²={k_squared}'
    )


def helmholtz_sin_k5():
    """
    H2: Medium frequency sine (-u'' + 25u = f)

    Exact: u(x) = sin(5πx)
    k² = 25 (k = 5)

    Medium frequency oscillations. Spectral methods should excel.
    """
    k = 5.0
    k_squared = k**2

    def exact(x):
        return np.sin(5 * np.pi * x)

    def forcing(x):
        # -u'' + k²*u = f
        return ((5 * np.pi)**2 + k_squared) * np.sin(5 * np.pi * x)

    return create_helmholtz_problem(
        name='Helmholtz_Sin_k5',
        domain=(0.0, 1.0),
        k_squared=k_squared,
        forcing=forcing,
        exact=exact,
        alpha=0.0,
        beta=0.0,
        description=f'Medium frequency sine, k²={k_squared}'
    )


def helmholtz_sin_k10():
    """
    H3: High frequency sine (-u'' + 100u = f)

    Exact: u(x) = sin(10πx)
    k² = 100 (k = 10)

    High frequency oscillations. Resolution requirements increase.
    """
    k = 10.0
    k_squared = k**2

    def exact(x):
        return np.sin(10 * np.pi * x)

    def forcing(x):
        # -u'' + k²*u = f
        return ((10 * np.pi)**2 + k_squared) * np.sin(10 * np.pi * x)

    return create_helmholtz_problem(
        name='Helmholtz_Sin_k10',
        domain=(0.0, 1.0),
        k_squared=k_squared,
        forcing=forcing,
        exact=exact,
        alpha=0.0,
        beta=0.0,
        description=f'High frequency sine, k²={k_squared}'
    )


def helmholtz_exp_k4():
    """
    H4: Exponential solution (-u'' + 16u = f)

    Exact: u(x) = e^(-2x) * sin(2πx)
    k² = 16 (k = 4)

    Tests rational approximation of exponential decay with oscillations.
    Rationals should naturally represent this well.
    """
    k = 4.0
    k_squared = k**2

    def exact(x):
        return np.exp(-2 * x) * np.sin(2 * np.pi * x)

    def forcing(x):
        # u = e^(-2x) * sin(2πx)
        # u' = -2*e^(-2x)*sin(2πx) + e^(-2x)*2π*cos(2πx)
        #    = e^(-2x)*[-2*sin(2πx) + 2π*cos(2πx)]
        # u'' = -2*e^(-2x)*[-2*sin(2πx) + 2π*cos(2πx)] +
        #       e^(-2x)*[-2*2π*cos(2πx) - 2π*2π*sin(2πx)]
        #     = e^(-2x)*[4*sin(2πx) - 4π*cos(2πx) - 4π*cos(2πx) - 4π²*sin(2πx)]
        #     = e^(-2x)*[(4 - 4π²)*sin(2πx) - 8π*cos(2πx)]
        u = exact(x)
        uxx = np.exp(-2 * x) * ((4 - 4 * np.pi**2) * np.sin(2 * np.pi * x) -
                                 8 * np.pi * np.cos(2 * np.pi * x))
        # -u'' + k²*u = f
        return -uxx + k_squared * u

    return create_helmholtz_problem(
        name='Helmholtz_Exp_k4',
        domain=(0.0, 1.0),
        k_squared=k_squared,
        forcing=forcing,
        exact=exact,
        alpha=exact(0.0),
        beta=exact(1.0),
        description=f'Exponential-oscillatory, k²={k_squared}, tests rational natural representation'
    )


def helmholtz_polynomial():
    """
    H5: Polynomial solution (-u'' + 4u = f)

    Exact: u(x) = x²(1-x)²
    k² = 4 (k = 2)

    Pure polynomial solution. Rationals offer no advantage; tests overhead.
    High-order FD and spectral should be extremely efficient.
    """
    k = 2.0
    k_squared = k**2

    def exact(x):
        return x**2 * (1 - x)**2

    def forcing(x):
        # u = x²(1-x)² = x²(1 - 2x + x²) = x² - 2x³ + x⁴
        # u' = 2x - 6x² + 4x³
        # u'' = 2 - 12x + 12x²
        uxx = 2 - 12 * x + 12 * x**2
        return -uxx + k_squared * exact(x)

    return create_helmholtz_problem(
        name='Helmholtz_Polynomial',
        domain=(0.0, 1.0),
        k_squared=k_squared,
        forcing=forcing,
        exact=exact,
        alpha=0.0,
        beta=0.0,
        description=f'Polynomial solution, k²={k_squared}, rationals have no advantage'
    )


def get_all_helmholtz_problems():
    """Return all Helmholtz test problems."""
    return [
        helmholtz_sin_k1(),
        helmholtz_sin_k5(),
        helmholtz_sin_k10(),
        helmholtz_exp_k4(),
        helmholtz_polynomial(),
    ]


if __name__ == "__main__":
    # Test problem definitions
    problems = get_all_helmholtz_problems()

    print("Helmholtz Test Problems")
    print("=" * 70)

    for prob in problems:
        print(f"\n{prob.name}:")
        print(f"  Domain: [{prob.domain[0]}, {prob.domain[1]}]")
        print(f"  k²: {prob.coeffs['k_squared']}")
        print(f"  BCs: u({prob.domain[0]}) = {prob.boundary_conditions['alpha']}, "
              f"u({prob.domain[1]}) = {prob.boundary_conditions['beta']}")
        print(f"  Description: {prob.description}")
        print(f"  Difficulty: {prob.expected_difficulty}")

        # Test that exact solution satisfies BCs
        u_a = prob.exact(prob.domain[0])
        u_b = prob.exact(prob.domain[1])
        print(f"  Verification: u(a) = {u_a:.6f} (should be {prob.boundary_conditions['alpha']:.6f})")
        print(f"                u(b) = {u_b:.6f} (should be {prob.boundary_conditions['beta']:.6f})")
