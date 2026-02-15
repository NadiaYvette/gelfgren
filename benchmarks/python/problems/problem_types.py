#!/usr/bin/env python3
"""
BVP Problem Type Definitions

Defines data structures for various types of boundary value problems.
"""

from dataclasses import dataclass
from typing import Callable, Tuple, Dict, Optional


@dataclass
class BVPProblem:
    """
    General boundary value problem specification.

    Represents BVPs of the form:
        L[u] = f(x)  on [a, b]
        u(a) = alpha, u(b) = beta

    Where L is a differential operator that can be:
    - Poisson: L[u] = -u''
    - Helmholtz: L[u] = -u'' + k²*u
    - Advection-Diffusion: L[u] = -ε*u'' + b*u'
    - Reaction-Diffusion: L[u] = -u'' + c*u
    - Variable Coefficient: L[u] = -(a(x)*u')'
    - General: L[u] = -ε*u'' + b*u' + c*u + k²*u
    """
    name: str
    """Problem name/identifier"""

    operator_type: str
    """Type of differential operator: 'poisson', 'helmholtz',
    'advection_diffusion', 'reaction_diffusion', 'variable_coefficient', 'general'"""

    domain: Tuple[float, float]
    """Domain boundaries (a, b)"""

    boundary_conditions: Dict[str, float]
    """Boundary values {'alpha': value at a, 'beta': value at b}"""

    forcing: Callable[[float], float]
    """Right-hand side forcing function f(x)"""

    exact: Optional[Callable[[float], float]]
    """Exact solution u(x) if known, None otherwise"""

    coeffs: Dict[str, any]
    """Operator-specific coefficients.

    For different operator types:
    - Poisson: {} (no extra coefficients)
    - Helmholtz: {'k_squared': k²}
    - Advection-Diffusion: {'epsilon': ε, 'b_coeff': b}
    - Reaction-Diffusion: {'c_coeff': c}
    - Variable Coefficient: {'a_func': a(x)}
    - General: {'epsilon': ε, 'b_coeff': b, 'c_coeff': c, 'k_squared': k²}
    """

    description: str = ""
    """Optional description of the problem"""

    expected_difficulty: str = "smooth"
    """Expected difficulty: 'smooth', 'boundary_layer', 'discontinuous',
    'oscillatory', 'singular'"""


def create_poisson_problem(name: str,
                          domain: Tuple[float, float],
                          forcing: Callable,
                          exact: Callable,
                          alpha: float = 0.0,
                          beta: float = 0.0,
                          description: str = "") -> BVPProblem:
    """
    Create a Poisson problem: -u'' = f(x).

    Parameters
    ----------
    name : str
        Problem identifier
    domain : tuple
        Domain boundaries (a, b)
    forcing : callable
        Right-hand side f(x)
    exact : callable
        Exact solution u(x)
    alpha, beta : float
        Boundary values
    description : str
        Problem description

    Returns
    -------
    BVPProblem
    """
    return BVPProblem(
        name=name,
        operator_type='poisson',
        domain=domain,
        boundary_conditions={'alpha': alpha, 'beta': beta},
        forcing=forcing,
        exact=exact,
        coeffs={},
        description=description,
        expected_difficulty='smooth'
    )


def create_helmholtz_problem(name: str,
                            domain: Tuple[float, float],
                            k_squared: float,
                            forcing: Callable,
                            exact: Callable,
                            alpha: float = 0.0,
                            beta: float = 0.0,
                            description: str = "") -> BVPProblem:
    """
    Create a Helmholtz problem: -u'' + k²*u = f(x).

    Parameters
    ----------
    name : str
        Problem identifier
    domain : tuple
        Domain boundaries (a, b)
    k_squared : float
        Helmholtz parameter k²
    forcing : callable
        Right-hand side f(x)
    exact : callable
        Exact solution u(x)
    alpha, beta : float
        Boundary values
    description : str
        Problem description

    Returns
    -------
    BVPProblem
    """
    difficulty = 'oscillatory' if k_squared > 50 else 'smooth'

    return BVPProblem(
        name=name,
        operator_type='helmholtz',
        domain=domain,
        boundary_conditions={'alpha': alpha, 'beta': beta},
        forcing=forcing,
        exact=exact,
        coeffs={'k_squared': k_squared},
        description=description,
        expected_difficulty=difficulty
    )


def create_advection_diffusion_problem(name: str,
                                      domain: Tuple[float, float],
                                      epsilon: float,
                                      b_coeff: float,
                                      forcing: Callable,
                                      exact: Callable,
                                      alpha: float = 0.0,
                                      beta: float = 0.0,
                                      description: str = "") -> BVPProblem:
    """
    Create an advection-diffusion problem: -ε*u'' + b*u' = f(x).

    Parameters
    ----------
    name : str
        Problem identifier
    domain : tuple
        Domain boundaries (a, b)
    epsilon : float
        Diffusion coefficient
    b_coeff : float
        Advection coefficient
    forcing : callable
        Right-hand side f(x)
    exact : callable
        Exact solution u(x)
    alpha, beta : float
        Boundary values
    description : str
        Problem description

    Returns
    -------
    BVPProblem
    """
    # Peclet number determines difficulty
    L = domain[1] - domain[0]
    Pe = abs(b_coeff) * L / epsilon if epsilon > 0 else float('inf')

    if Pe > 100:
        difficulty = 'boundary_layer'
    elif Pe > 10:
        difficulty = 'sharp_gradient'
    else:
        difficulty = 'smooth'

    return BVPProblem(
        name=name,
        operator_type='advection_diffusion',
        domain=domain,
        boundary_conditions={'alpha': alpha, 'beta': beta},
        forcing=forcing,
        exact=exact,
        coeffs={'epsilon': epsilon, 'b_coeff': b_coeff},
        description=description,
        expected_difficulty=difficulty
    )


def create_reaction_diffusion_problem(name: str,
                                     domain: Tuple[float, float],
                                     c_coeff: float,
                                     forcing: Callable,
                                     exact: Callable,
                                     alpha: float = 0.0,
                                     beta: float = 0.0,
                                     description: str = "") -> BVPProblem:
    """
    Create a reaction-diffusion problem: -u'' + c*u = f(x).

    Parameters
    ----------
    name : str
        Problem identifier
    domain : tuple
        Domain boundaries (a, b)
    c_coeff : float
        Reaction coefficient
    forcing : callable
        Right-hand side f(x)
    exact : callable
        Exact solution u(x)
    alpha, beta : float
        Boundary values
    description : str
        Problem description

    Returns
    -------
    BVPProblem
    """
    difficulty = 'smooth' if c_coeff < 50 else 'boundary_layer'

    return BVPProblem(
        name=name,
        operator_type='reaction_diffusion',
        domain=domain,
        boundary_conditions={'alpha': alpha, 'beta': beta},
        forcing=forcing,
        exact=exact,
        coeffs={'c_coeff': c_coeff},
        description=description,
        expected_difficulty=difficulty
    )


def create_variable_coefficient_problem(name: str,
                                       domain: Tuple[float, float],
                                       a_func: Callable,
                                       forcing: Callable,
                                       exact: Callable,
                                       alpha: float = 0.0,
                                       beta: float = 0.0,
                                       description: str = "",
                                       expected_difficulty: str = "smooth") -> BVPProblem:
    """
    Create a variable coefficient problem: -(a(x)*u')' = f(x).

    Parameters
    ----------
    name : str
        Problem identifier
    domain : tuple
        Domain boundaries (a, b)
    a_func : callable
        Variable diffusion coefficient a(x)
    forcing : callable
        Right-hand side f(x)
    exact : callable
        Exact solution u(x)
    alpha, beta : float
        Boundary values
    description : str
        Problem description
    expected_difficulty : str
        Expected difficulty level

    Returns
    -------
    BVPProblem
    """
    return BVPProblem(
        name=name,
        operator_type='variable_coefficient',
        domain=domain,
        boundary_conditions={'alpha': alpha, 'beta': beta},
        forcing=forcing,
        exact=exact,
        coeffs={'a_func': a_func},
        description=description,
        expected_difficulty=expected_difficulty
    )
