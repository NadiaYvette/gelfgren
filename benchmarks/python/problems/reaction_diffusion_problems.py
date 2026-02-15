#!/usr/bin/env python3
"""
Reaction-Diffusion BVP Test Problems

Test problems for: -u'' + c*u = f(x) on [0, 1]
                   u(0) = α, u(1) = β

These problems test method performance on exponential solutions.
Rationals should excel when c > 0 (natural exponential representation).
"""

import numpy as np
from problem_types import create_reaction_diffusion_problem


def reaction_diffusion_mild():
    """
    RD1: Mild reaction (c = 1)

    Exact: u(x) = sin(πx)
    c = 1

    Smooth oscillatory solution with mild reaction term.
    All spectral methods should converge exponentially.
    """
    c = 1.0

    def exact(x):
        return np.sin(np.pi * x)

    def forcing(x):
        # -u'' + c*u = f
        # -(-π²*sin(πx)) + c*sin(πx) = f
        return (np.pi**2 + c) * np.sin(np.pi * x)

    return create_reaction_diffusion_problem(
        name='ReacDiff_Mild_c1',
        domain=(0.0, 1.0),
        c_coeff=c,
        forcing=forcing,
        exact=exact,
        alpha=0.0,
        beta=0.0,
        description=f'Mild reaction, c={c}, smooth oscillatory'
    )


def reaction_diffusion_exponential_layers():
    """
    RD2: Exponential boundary layers (c = 10)

    Exact: u(x) = sinh(√c * x) / sinh(√c)
    c = 10, √c ≈ 3.16

    Solution: -u'' + c*u = 0 with u(0)=0, u(1)=1
    Creates exponential boundary layers at both endpoints.
    Rationals should represent this naturally well.
    """
    c = 10.0
    sqrt_c = np.sqrt(c)

    def exact(x):
        return np.sinh(sqrt_c * x) / np.sinh(sqrt_c)

    def forcing(x):
        # -u'' + c*u = 0 for this solution
        return np.zeros_like(x) if hasattr(x, '__len__') else 0.0

    return create_reaction_diffusion_problem(
        name='ReacDiff_Exp_c10',
        domain=(0.0, 1.0),
        c_coeff=c,
        forcing=forcing,
        exact=exact,
        alpha=0.0,
        beta=1.0,
        description=f'Exponential layers, c={c}, rationals should excel'
    )


def reaction_diffusion_steep():
    """
    RD3: Steep exponential (c = 100)

    Exact: u(x) = sinh(10x) / sinh(10)
    c = 100, √c = 10

    Very steep exponential growth with boundary layers.
    Rationals should still represent naturally, but resolution demands high.
    This is smooth (unlike discontinuous) - rationals have advantage.
    """
    c = 100.0
    sqrt_c = np.sqrt(c)

    def exact(x):
        return np.sinh(sqrt_c * x) / np.sinh(sqrt_c)

    def forcing(x):
        return np.zeros_like(x) if hasattr(x, '__len__') else 0.0

    return create_reaction_diffusion_problem(
        name='ReacDiff_Steep_c100',
        domain=(0.0, 1.0),
        c_coeff=c,
        forcing=forcing,
        exact=exact,
        alpha=0.0,
        beta=1.0,
        description=f'Steep exponential, c={c}, smooth but demanding, rationals should excel'
    )


def reaction_diffusion_oscillatory():
    """
    RD4: Variable reaction coefficient c(x)

    Exact: u(x) = e^(-x) * cos(5πx)
    c(x) = variable to create forcing

    Actually, we'll treat this as constant coefficient with
    special forcing that creates oscillatory-exponential solution.
    Tests rational representation of damped oscillations.
    """
    c = 5.0

    def exact(x):
        return np.exp(-x) * np.cos(5 * np.pi * x)

    def forcing(x):
        # u = e^(-x) * cos(5πx)
        # u' = -e^(-x)*cos(5πx) - 5π*e^(-x)*sin(5πx)
        #    = e^(-x)*[-cos(5πx) - 5π*sin(5πx)]
        # u'' = -e^(-x)*[-cos(5πx) - 5π*sin(5πx)] +
        #       e^(-x)*[5π*sin(5πx) - 25π²*cos(5πx)]
        #     = e^(-x)*[cos(5πx) + 5π*sin(5πx) + 5π*sin(5πx) - 25π²*cos(5πx)]
        #     = e^(-x)*[(1 - 25π²)*cos(5πx) + 10π*sin(5πx)]

        u = exact(x)
        uxx = np.exp(-x) * ((1 - 25 * np.pi**2) * np.cos(5 * np.pi * x) +
                             10 * np.pi * np.sin(5 * np.pi * x))
        return -uxx + c * u

    return create_reaction_diffusion_problem(
        name='ReacDiff_Oscillatory_c5',
        domain=(0.0, 1.0),
        c_coeff=c,
        forcing=forcing,
        exact=exact,
        alpha=exact(0.0),
        beta=exact(1.0),
        description=f'Damped oscillatory, c={c}, tests rational exponential representation'
    )


def reaction_diffusion_polynomial():
    """
    RD5: Polynomial solution (c = 4)

    Exact: u(x) = x²(1-x)² = x² - 2x³ + x⁴
    c = 4

    Pure polynomial. Rationals offer no advantage; tests overhead.
    High-order polynomial methods should be most efficient.
    """
    c = 4.0

    def exact(x):
        return x**2 * (1 - x)**2

    def forcing(x):
        # u = x² - 2x³ + x⁴
        # u' = 2x - 6x² + 4x³
        # u'' = 2 - 12x + 12x²
        uxx = 2 - 12 * x + 12 * x**2
        return -uxx + c * exact(x)

    return create_reaction_diffusion_problem(
        name='ReacDiff_Polynomial_c4',
        domain=(0.0, 1.0),
        c_coeff=c,
        forcing=forcing,
        exact=exact,
        alpha=0.0,
        beta=0.0,
        description=f'Polynomial solution, c={c}, rationals have no advantage'
    )


def get_all_reaction_diffusion_problems():
    """Return all reaction-diffusion test problems."""
    return [
        reaction_diffusion_mild(),
        reaction_diffusion_exponential_layers(),
        reaction_diffusion_steep(),
        reaction_diffusion_oscillatory(),
        reaction_diffusion_polynomial(),
    ]


if __name__ == "__main__":
    # Test problem definitions
    problems = get_all_reaction_diffusion_problems()

    print("Reaction-Diffusion Test Problems")
    print("=" * 70)

    for prob in problems:
        print(f"\n{prob.name}:")
        print(f"  Domain: [{prob.domain[0]}, {prob.domain[1]}]")
        print(f"  c: {prob.coeffs['c_coeff']}")
        print(f"  BCs: u({prob.domain[0]}) = {prob.boundary_conditions['alpha']:.6f}, "
              f"u({prob.domain[1]}) = {prob.boundary_conditions['beta']:.6f}")
        print(f"  Description: {prob.description}")
        print(f"  Difficulty: {prob.expected_difficulty}")

        # Test that exact solution satisfies BCs
        u_a = prob.exact(prob.domain[0])
        u_b = prob.exact(prob.domain[1])
        print(f"  Verification: u(a) = {u_a:.6f} (should be {prob.boundary_conditions['alpha']:.6f})")
        print(f"                u(b) = {u_b:.6f} (should be {prob.boundary_conditions['beta']:.6f})")
