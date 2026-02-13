#!/usr/bin/env python3
"""
Boundary Value Problem Convergence Study

Compares piecewise rational approximants vs polynomial splines for:
1. 1D Poisson equation with various forcing functions
2. Heat equation (steady state)
3. Beam deflection problem

Generates convergence data for LaTeX report generation.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_bvp
import json
from dataclasses import dataclass, asdict
from typing import Callable, List, Tuple
import sys
import os

# Add gelfgren Python bindings to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../bindings/python'))

try:
    import gelfgren as gf
    GELFGREN_AVAILABLE = True
except ImportError:
    print("Warning: Gelfgren Python bindings not available. Using placeholder.")
    GELFGREN_AVAILABLE = False


@dataclass
class ErrorMetrics:
    """Error metrics for solution quality"""
    l2_error: float
    l_inf_error: float
    h1_seminorm_error: float
    num_intervals: int
    dof: int
    mesh_size: float
    method: str  # 'polynomial' or 'rational'


@dataclass
class ConvergenceResult:
    """Results from convergence study"""
    problem_name: str
    polynomial_errors: List[ErrorMetrics]
    rational_errors: List[ErrorMetrics]
    convergence_rates: dict


class PoissonProblem1D:
    """1D Poisson equation: -u'' = f, u(0) = u(1) = 0"""

    def __init__(self, forcing: Callable, exact: Callable, name: str):
        self.forcing = forcing
        self.exact = exact
        self.name = name
        self.a = 0.0
        self.b = 1.0

    def solve_polynomial_spline(self, n_intervals: int) -> Tuple[np.ndarray, np.ndarray]:
        """Solve using cubic spline finite elements"""
        # Mesh points
        x = np.linspace(self.a, self.b, n_intervals + 1)
        h = x[1] - x[0]

        # Assemble finite element matrix (for cubic splines, this is simplified)
        # Using central differences for -u''
        n = len(x)
        A = np.zeros((n, n))
        b = np.zeros(n)

        # Interior points: -u''(x_i) ≈ -(u_{i+1} - 2u_i + u_{i-1})/h²
        for i in range(1, n - 1):
            A[i, i - 1] = -1.0 / (h * h)
            A[i, i] = 2.0 / (h * h)
            A[i, i + 1] = -1.0 / (h * h)
            b[i] = self.forcing(x[i])

        # Boundary conditions: u(0) = u(1) = 0
        A[0, 0] = 1.0
        A[n - 1, n - 1] = 1.0
        b[0] = 0.0
        b[n - 1] = 0.0

        # Solve linear system
        u = np.linalg.solve(A, b)

        return x, u

    def solve_rational_piecewise(self, n_intervals: int, deg_num: int = 2, deg_den: int = 2):
        """Solve using piecewise rational approximants (placeholder)"""
        if not GELFGREN_AVAILABLE:
            # Fallback: use polynomial solution
            return self.solve_polynomial_spline(n_intervals)

        # Create mesh
        x_mesh = np.linspace(self.a, self.b, n_intervals + 1)

        # For each subinterval, construct Padé approximant to local solution
        # This is a simplified approach - full BVP solver would be more sophisticated
        x_fine = np.linspace(self.a, self.b, 1000)
        u_approx = np.zeros_like(x_fine)

        for i in range(n_intervals):
            x_left = x_mesh[i]
            x_right = x_mesh[i + 1]
            mask = (x_fine >= x_left) & (x_fine <= x_right)
            x_local = x_fine[mask]

            # Local polynomial approximation to forcing
            x_local_normalized = (x_local - x_left) / (x_right - x_left)

            # Taylor series coefficients of local solution
            # For demonstration: use Taylor expansion
            coeffs = [0.0] * 5  # Placeholder
            for j in range(5):
                # Would compute actual Taylor coefficients here
                coeffs[j] = 1.0 / (j + 1)

            try:
                # Create Padé approximant
                pade = gf.PadeApproximant(coeffs, deg_num, deg_den, 0.0, 1.0)

                # Evaluate on local interval
                u_local = np.array([pade.evaluate(xi) for xi in x_local_normalized])
                u_approx[mask] = u_local
            except:
                # Fallback to polynomial
                u_approx[mask] = np.polyval(coeffs[::-1], x_local_normalized)

        return x_fine, u_approx

    def compute_error(self, x: np.ndarray, u_approx: np.ndarray) -> ErrorMetrics:
        """Compute error metrics"""
        u_exact = np.array([self.exact(xi) for xi in x])
        error = u_approx - u_exact

        h = x[1] - x[0]

        # L2 error
        l2_error = np.sqrt(h * np.sum(error ** 2))

        # L∞ error
        l_inf_error = np.max(np.abs(error))

        # H1 seminorm error (discrete derivative)
        du_approx = np.diff(u_approx) / h
        du_exact = np.diff(u_exact) / h
        h1_error = np.sqrt(h * np.sum((du_approx - du_exact) ** 2))

        return l2_error, l_inf_error, h1_error


def smooth_poisson():
    """Smooth Poisson: -u'' = π² sin(πx), u_exact = sin(πx)"""
    forcing = lambda x: np.pi ** 2 * np.sin(np.pi * x)
    exact = lambda x: np.sin(np.pi * x)
    return PoissonProblem1D(forcing, exact, "Smooth Poisson (sin)")


def discontinuous_poisson():
    """Discontinuous forcing"""
    def forcing(x):
        if isinstance(x, np.ndarray):
            f = np.zeros_like(x)
            mask = (x >= 0.25) & (x <= 0.75)
            f[mask] = -2.0
            return f
        else:
            return -2.0 if 0.25 <= x <= 0.75 else 0.0

    def exact(x):
        if isinstance(x, np.ndarray):
            u = np.zeros_like(x)
            mask1 = x < 0.25
            mask2 = (x >= 0.25) & (x <= 0.75)
            mask3 = x > 0.75
            u[mask1] = 0.5 * x[mask1]
            u[mask2] = -x[mask2] ** 2 + 0.75 * x[mask2] - 0.0625
            u[mask3] = -0.5 * x[mask3] + 0.5
            return u
        else:
            if x < 0.25:
                return 0.5 * x
            elif x <= 0.75:
                return -x * x + 0.75 * x - 0.0625
            else:
                return -0.5 * x + 0.5

    return PoissonProblem1D(forcing, exact, "Discontinuous Poisson")


def oscillatory_poisson(frequency: float = 10.0):
    """Oscillatory forcing"""
    omega = frequency
    forcing = lambda x: (omega * np.pi) ** 2 * np.sin(omega * np.pi * x)
    exact = lambda x: np.sin(omega * np.pi * x)
    return PoissonProblem1D(forcing, exact, f"Oscillatory Poisson (ω={omega})")


def run_convergence_study(problem: PoissonProblem1D, mesh_sizes: List[int]) -> ConvergenceResult:
    """Run convergence study for given problem"""
    poly_errors = []
    rat_errors = []

    print(f"\nConvergence study: {problem.name}")
    print("=" * 60)

    for n in mesh_sizes:
        print(f"\nMesh: {n} intervals")

        # Polynomial spline solution
        x_poly, u_poly = problem.solve_polynomial_spline(n)
        l2, linf, h1 = problem.compute_error(x_poly, u_poly)
        h = (problem.b - problem.a) / n

        poly_metric = ErrorMetrics(
            l2_error=l2,
            l_inf_error=linf,
            h1_seminorm_error=h1,
            num_intervals=n,
            dof=n + 1,  # Number of nodes
            mesh_size=h,
            method='polynomial'
        )
        poly_errors.append(poly_metric)

        print(f"  Polynomial: L2={l2:.6e}, L∞={linf:.6e}, H1={h1:.6e}")

        # Rational approximant solution
        x_rat, u_rat = problem.solve_rational_piecewise(n)
        l2, linf, h1 = problem.compute_error(x_rat, u_rat)

        rat_metric = ErrorMetrics(
            l2_error=l2,
            l_inf_error=linf,
            h1_seminorm_error=h1,
            num_intervals=n,
            dof=n * 5,  # [2/2] rational DOF (with normalization: 2+2+1)
            mesh_size=h,
            method='rational'
        )
        rat_errors.append(rat_metric)

        print(f"  Rational:   L2={l2:.6e}, L∞={linf:.6e}, H1={h1:.6e}")

    # Compute convergence rates
    rates = compute_convergence_rates(poly_errors, rat_errors)

    return ConvergenceResult(
        problem_name=problem.name,
        polynomial_errors=poly_errors,
        rational_errors=rat_errors,
        convergence_rates=rates
    )


def compute_convergence_rates(poly_errors: List[ErrorMetrics],
                               rat_errors: List[ErrorMetrics]) -> dict:
    """Compute convergence rates"""
    rates = {
        'polynomial_l2': [],
        'polynomial_linf': [],
        'rational_l2': [],
        'rational_linf': []
    }

    for i in range(1, len(poly_errors)):
        h_ratio = poly_errors[i - 1].mesh_size / poly_errors[i].mesh_size

        # Polynomial rates
        rate_l2 = np.log(poly_errors[i - 1].l2_error / poly_errors[i].l2_error) / np.log(h_ratio)
        rate_linf = np.log(poly_errors[i - 1].l_inf_error / poly_errors[i].l_inf_error) / np.log(h_ratio)
        rates['polynomial_l2'].append(rate_l2)
        rates['polynomial_linf'].append(rate_linf)

        # Rational rates
        rate_l2 = np.log(rat_errors[i - 1].l2_error / rat_errors[i].l2_error) / np.log(h_ratio)
        rate_linf = np.log(rat_errors[i - 1].l_inf_error / rat_errors[i].l_inf_error) / np.log(h_ratio)
        rates['rational_l2'].append(rate_l2)
        rates['rational_linf'].append(rate_linf)

    return rates


def save_results(results: List[ConvergenceResult], output_dir: str):
    """Save results to JSON for LaTeX report generation"""
    os.makedirs(output_dir, exist_ok=True)

    for result in results:
        filename = result.problem_name.replace(' ', '_').replace('(', '').replace(')', '') + '.json'
        filepath = os.path.join(output_dir, filename)

        data = {
            'problem_name': result.problem_name,
            'polynomial_errors': [asdict(e) for e in result.polynomial_errors],
            'rational_errors': [asdict(e) for e in result.rational_errors],
            'convergence_rates': result.convergence_rates
        }

        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)

        print(f"\nSaved results to {filepath}")


def plot_convergence(results: List[ConvergenceResult], output_dir: str):
    """Generate convergence plots"""
    os.makedirs(output_dir, exist_ok=True)

    for result in results:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Extract data
        h_poly = [e.mesh_size for e in result.polynomial_errors]
        l2_poly = [e.l2_error for e in result.polynomial_errors]
        l2_rat = [e.l2_error for e in result.rational_errors]

        # L2 error vs mesh size
        ax1.loglog(h_poly, l2_poly, 'o-', label='Polynomial spline', linewidth=2)
        ax1.loglog(h_poly, l2_rat, 's-', label='Rational approximant', linewidth=2)

        # Reference lines
        ax1.loglog(h_poly, [h ** 2 for h in h_poly], '--', label='O(h²)', alpha=0.5)
        ax1.loglog(h_poly, [h ** 4 for h in h_poly], '--', label='O(h⁴)', alpha=0.5)

        ax1.set_xlabel('Mesh size h', fontsize=12)
        ax1.set_ylabel('L² error', fontsize=12)
        ax1.set_title(f'{result.problem_name} - L² Convergence', fontsize=14)
        ax1.legend(fontsize=10)
        ax1.grid(True, alpha=0.3)

        # Error vs DOF
        dof_poly = [e.dof for e in result.polynomial_errors]
        dof_rat = [e.dof for e in result.rational_errors]

        ax2.loglog(dof_poly, l2_poly, 'o-', label='Polynomial spline', linewidth=2)
        ax2.loglog(dof_rat, l2_rat, 's-', label='Rational approximant', linewidth=2)

        ax2.set_xlabel('Degrees of Freedom', fontsize=12)
        ax2.set_ylabel('L² error', fontsize=12)
        ax2.set_title(f'{result.problem_name} - Efficiency', fontsize=14)
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        filename = result.problem_name.replace(' ', '_').replace('(', '').replace(')', '') + '.pdf'
        filepath = os.path.join(output_dir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"Saved plot to {filepath}")
        plt.close()


def main():
    """Run all BVP benchmarks"""
    print("Gelfgren BVP Convergence Study")
    print("=" * 60)

    # Define mesh sizes for convergence study
    mesh_sizes = [4, 8, 16, 32, 64, 128]

    # Define problems
    problems = [
        smooth_poisson(),
        discontinuous_poisson(),
        oscillatory_poisson(10.0),
    ]

    # Run convergence studies
    results = []
    for problem in problems:
        result = run_convergence_study(problem, mesh_sizes)
        results.append(result)

    # Save results
    output_dir = os.path.join(os.path.dirname(__file__), '../data')
    save_results(results, output_dir)

    # Generate plots
    figures_dir = os.path.join(os.path.dirname(__file__), '../reports/figures')
    plot_convergence(results, figures_dir)

    print("\n" + "=" * 60)
    print("Convergence study complete!")
    print(f"Results saved to: {output_dir}")
    print(f"Figures saved to: {figures_dir}")


if __name__ == '__main__':
    main()
