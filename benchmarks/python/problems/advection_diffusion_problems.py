#!/usr/bin/env python3
"""
Advection-Diffusion BVP Test Problems

Test problems for: -ε*u'' + b*u' = f(x) on [0, 1]
                   u(0) = α, u(1) = β

These problems test method performance on boundary layers.
CRITICAL: Rationals expected to catastrophically fail at ε ≤ 0.01.
"""

import numpy as np
from problem_types import create_advection_diffusion_problem


def advection_diffusion_mild():
    """
    AD1: Mild boundary layer (ε = 0.1, b = 1)

    Exact: u(x) = (exp(-x/ε) - exp(-1/ε)) / (1 - exp(-1/ε))
    Peclet number: Pe = b*L/ε = 1*1/0.1 = 10

    Boundary layer thickness ~ ε = 0.1
    Smooth enough that rationals should work, though not optimally.
    """
    epsilon = 0.1
    b = 1.0

    def exact(x):
        return (np.exp(-x / epsilon) - np.exp(-1.0 / epsilon)) / (1 - np.exp(-1.0 / epsilon))

    def forcing(x):
        # Since -ε*u'' + b*u' = 0 for this solution, f = 0
        return np.zeros_like(x) if hasattr(x, '__len__') else 0.0

    return create_advection_diffusion_problem(
        name='AdvDiff_Mild_eps0.1',
        domain=(0.0, 1.0),
        epsilon=epsilon,
        b_coeff=b,
        forcing=forcing,
        exact=exact,
        alpha=1.0,
        beta=0.0,
        description=f'Mild boundary layer, ε={epsilon}, Pe=10, rationals may struggle slightly'
    )


def advection_diffusion_noticeable():
    """
    AD2: Noticeable boundary layer (ε = 0.01, b = 1)

    Exact: u(x) = (exp(-x/ε) - exp(-1/ε)) / (1 - exp(-1/ε))
    Peclet number: Pe = 100

    Boundary layer thickness ~ ε = 0.01
    Rationals should start to degrade significantly.
    """
    epsilon = 0.01
    b = 1.0

    def exact(x):
        return (np.exp(-x / epsilon) - np.exp(-1.0 / epsilon)) / (1 - np.exp(-1.0 / epsilon))

    def forcing(x):
        return np.zeros_like(x) if hasattr(x, '__len__') else 0.0

    return create_advection_diffusion_problem(
        name='AdvDiff_Noticeable_eps0.01',
        domain=(0.0, 1.0),
        epsilon=epsilon,
        b_coeff=b,
        forcing=forcing,
        exact=exact,
        alpha=1.0,
        beta=0.0,
        description=f'Noticeable boundary layer, ε={epsilon}, Pe=100, rationals degrade'
    )


def advection_diffusion_sharp():
    """
    AD3: Sharp boundary layer (ε = 0.001, b = 1)

    Exact: u(x) = (exp(-x/ε) - exp(-1/ε)) / (1 - exp(-1/ε))
    Peclet number: Pe = 1000

    Boundary layer thickness ~ ε = 0.001
    **CRITICAL: Rationals expected to CATASTROPHICALLY FAIL**

    This is a key test case - rationals cannot represent sharp gradients
    without poles, leading to complete breakdown.
    """
    epsilon = 0.001
    b = 1.0

    def exact(x):
        return (np.exp(-x / epsilon) - np.exp(-1.0 / epsilon)) / (1 - np.exp(-1.0 / epsilon))

    def forcing(x):
        return np.zeros_like(x) if hasattr(x, '__len__') else 0.0

    return create_advection_diffusion_problem(
        name='AdvDiff_Sharp_eps0.001',
        domain=(0.0, 1.0),
        epsilon=epsilon,
        b_coeff=b,
        forcing=forcing,
        exact=exact,
        alpha=1.0,
        beta=0.0,
        description=f'Sharp boundary layer, ε={epsilon}, Pe=1000, **RATIONALS FAIL CATASTROPHICALLY**'
    )


def advection_diffusion_interior_layer():
    """
    AD4: Interior layer (turning point problem)

    Exact: u(x) = x + εln((1+exp(-(x-0.5)/ε))/(1+exp(0.5/ε)))
    Advection coefficient: b(x) = 2x - 1 (changes sign at x=0.5)
    ε = 0.01

    Interior layer at x = 0.5 where advection changes direction.
    Extremely challenging for all methods without adaptive refinement.
    """
    epsilon = 0.01

    def b_func(x):
        return 2 * x - 1

    # For this problem, we use a simplified manufactured solution
    # that creates an interior layer
    def exact(x):
        # Smooth transition with sharp gradient at x=0.5
        return 0.5 * (1 + np.tanh((x - 0.5) / (2 * epsilon)))

    def forcing(x):
        # Compute f from -ε*u'' + b(x)*u' = f
        # u = 0.5*(1 + tanh((x-0.5)/(2ε)))
        # u' = 0.5 * sech²((x-0.5)/(2ε)) / (2ε) = sech²(s)/(4ε) where s=(x-0.5)/(2ε)
        # u'' = -2*sech²(s)*tanh(s)/(4ε)² = -sech²(s)*tanh(s)/(8ε³)

        s = (x - 0.5) / (2 * epsilon)
        sech = 1.0 / np.cosh(s)
        sech2 = sech**2
        tanh = np.tanh(s)

        ux = sech2 / (4 * epsilon)
        uxx = -sech2 * tanh / (8 * epsilon**3)

        b_val = b_func(x)
        return -epsilon * uxx + b_val * ux

    return create_advection_diffusion_problem(
        name='AdvDiff_Interior_eps0.01',
        domain=(0.0, 1.0),
        epsilon=epsilon,
        b_coeff=0.0,  # Variable, but we set to 0 and handle separately
        forcing=forcing,
        exact=exact,
        alpha=exact(0.0),
        beta=exact(1.0),
        description=f'Interior layer at x=0.5, ε={epsilon}, extremely challenging'
    )


def advection_diffusion_smooth_dominated():
    """
    AD5: Smooth advection-dominated (ε = 0.1, no layer)

    Exact: u(x) = x(1-x) (polynomial)
    ε = 0.1, chosen forcing to avoid boundary layer

    Advection-dominated but smooth throughout.
    All methods should work well; tests overhead of different approaches.
    """
    epsilon = 0.1
    b = 1.0

    def exact(x):
        return x * (1 - x)

    def forcing(x):
        # u = x(1-x) = x - x²
        # u' = 1 - 2x
        # u'' = -2
        # -ε*u'' + b*u' = -ε*(-2) + b*(1-2x) = 2ε + b - 2bx
        return 2 * epsilon + b * (1 - 2 * x)

    return create_advection_diffusion_problem(
        name='AdvDiff_Smooth_eps0.1',
        domain=(0.0, 1.0),
        epsilon=epsilon,
        b_coeff=b,
        forcing=forcing,
        exact=exact,
        alpha=0.0,
        beta=0.0,
        description=f'Smooth advection-dominated, ε={epsilon}, no boundary layer'
    )


def get_all_advection_diffusion_problems():
    """Return all advection-diffusion test problems."""
    return [
        advection_diffusion_mild(),
        advection_diffusion_noticeable(),
        advection_diffusion_sharp(),  # CRITICAL failure test
        advection_diffusion_interior_layer(),
        advection_diffusion_smooth_dominated(),
    ]


if __name__ == "__main__":
    # Test problem definitions
    problems = get_all_advection_diffusion_problems()

    print("Advection-Diffusion Test Problems")
    print("=" * 70)

    for prob in problems:
        print(f"\n{prob.name}:")
        print(f"  Domain: [{prob.domain[0]}, {prob.domain[1]}]")
        print(f"  ε: {prob.coeffs['epsilon']}, b: {prob.coeffs['b_coeff']}")
        if prob.coeffs['epsilon'] > 0:
            L = prob.domain[1] - prob.domain[0]
            Pe = abs(prob.coeffs['b_coeff']) * L / prob.coeffs['epsilon']
            print(f"  Peclet number: Pe = {Pe:.1f}")
        print(f"  BCs: u({prob.domain[0]}) = {prob.boundary_conditions['alpha']:.6f}, "
              f"u({prob.domain[1]}) = {prob.boundary_conditions['beta']:.6f}")
        print(f"  Description: {prob.description}")
        print(f"  Difficulty: {prob.expected_difficulty}")

        # Test that exact solution satisfies BCs
        u_a = prob.exact(prob.domain[0])
        u_b = prob.exact(prob.domain[1])
        print(f"  Verification: u(a) = {u_a:.6f} (should be {prob.boundary_conditions['alpha']:.6f})")
        print(f"                u(b) = {u_b:.6f} (should be {prob.boundary_conditions['beta']:.6f})")
