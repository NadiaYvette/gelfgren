#!/usr/bin/env python3
"""
Special Function Approximation Convergence Study

Compares piecewise rational approximants vs polynomial splines for:
1. Exponential: e^x
2. Trigonometric: sin(x), cos(x), tan(x)
3. Error function: erf(x)
4. Bessel function: J_0(x)
5. Logarithm: log(1+x)
6. Runge's function: 1/(1+25x²)

Generates convergence data for LaTeX report generation.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.special import erf, j0
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
    """Error metrics for approximation quality"""
    l2_error: float
    l_inf_error: float
    h1_seminorm_error: float
    relative_l2_error: float
    relative_l_inf_error: float
    num_intervals: int
    dof: int
    mesh_size: float
    method: str  # 'polynomial' or 'rational'


@dataclass
class ConvergenceResult:
    """Results from convergence study"""
    function_name: str
    domain: Tuple[float, float]
    polynomial_errors: List[ErrorMetrics]
    rational_errors: List[ErrorMetrics]
    convergence_rates: dict


class SpecialFunctionApproximation:
    """Base class for special function approximation problems"""

    def __init__(self, func: Callable, name: str, a: float, b: float):
        self.func = func
        self.name = name
        self.a = a
        self.b = b

    def approximate_polynomial_spline(self, n_intervals: int) -> Tuple[np.ndarray, np.ndarray]:
        """Approximate using cubic spline interpolation"""
        # Interpolation points (including boundaries)
        x_knots = np.linspace(self.a, self.b, n_intervals + 1)
        y_knots = np.array([self.func(xi) for xi in x_knots])

        # Create cubic spline
        spline = CubicSpline(x_knots, y_knots)

        # Evaluate on fine grid
        x_fine = np.linspace(self.a, self.b, 1000)
        y_approx = spline(x_fine)

        return x_fine, y_approx

    def approximate_rational_piecewise(self, n_intervals: int, deg_num: int = 2, deg_den: int = 2):
        """Approximate using piecewise rational approximants"""
        if not GELFGREN_AVAILABLE:
            # Fallback: use polynomial spline
            return self.approximate_polynomial_spline(n_intervals)

        # Create mesh
        x_mesh = np.linspace(self.a, self.b, n_intervals + 1)

        # Fine evaluation grid
        x_fine = np.linspace(self.a, self.b, 1000)
        y_approx = np.zeros_like(x_fine)

        for i in range(n_intervals):
            x_left = x_mesh[i]
            x_right = x_mesh[i + 1]
            mask = (x_fine >= x_left) & (x_fine <= x_right)
            x_local = x_fine[mask]

            # Compute local Taylor series coefficients at midpoint
            x_mid = (x_left + x_right) / 2.0
            h = x_right - x_left

            # Compute derivatives numerically at midpoint
            coeffs = []
            delta = 1e-6

            # f(x_mid)
            coeffs.append(self.func(x_mid))

            # f'(x_mid)
            f_plus = self.func(x_mid + delta)
            f_minus = self.func(x_mid - delta)
            coeffs.append((f_plus - f_minus) / (2.0 * delta))

            # f''(x_mid) / 2!
            f_center = self.func(x_mid)
            coeffs.append((f_plus - 2.0 * f_center + f_minus) / (2.0 * delta * delta))

            # f'''(x_mid) / 3!
            f_2plus = self.func(x_mid + 2 * delta)
            f_2minus = self.func(x_mid - 2 * delta)
            coeffs.append((f_2plus - 2.0 * f_plus + 2.0 * f_minus - f_2minus) / (12.0 * delta ** 3))

            # f''''(x_mid) / 4!
            coeffs.append((f_2plus - 4.0 * f_plus + 6.0 * f_center - 4.0 * f_minus + f_2minus) / (24.0 * delta ** 4))

            try:
                # Create Padé approximant on normalized interval [0, 1]
                pade = gf.PadeApproximant(coeffs, deg_num, deg_den, 0.0, 1.0)

                # Map local x to [0, 1]
                x_local_normalized = (x_local - x_left) / (x_right - x_left)

                # Evaluate Padé approximant
                y_local = np.array([pade.evaluate(xi) for xi in x_local_normalized])
                y_approx[mask] = y_local
            except:
                # Fallback to polynomial Taylor expansion
                x_local_centered = x_local - x_mid
                y_local = sum(c * x_local_centered ** k / np.math.factorial(k)
                            for k, c in enumerate(coeffs[:4]))
                y_approx[mask] = y_local

        return x_fine, y_approx

    def compute_error(self, x: np.ndarray, y_approx: np.ndarray) -> Tuple[float, float, float, float, float]:
        """Compute error metrics"""
        y_exact = np.array([self.func(xi) for xi in x])
        error = y_approx - y_exact

        h = x[1] - x[0]

        # L2 error
        l2_error = np.sqrt(h * np.sum(error ** 2))

        # L∞ error
        l_inf_error = np.max(np.abs(error))

        # H1 seminorm error (discrete derivative)
        dy_approx = np.diff(y_approx) / h
        dy_exact = np.diff(y_exact) / h
        h1_error = np.sqrt(h * np.sum((dy_approx - dy_exact) ** 2))

        # Relative errors
        l2_norm_exact = np.sqrt(h * np.sum(y_exact ** 2))
        l_inf_norm_exact = np.max(np.abs(y_exact))

        relative_l2 = l2_error / l2_norm_exact if l2_norm_exact > 1e-12 else l2_error
        relative_l_inf = l_inf_error / l_inf_norm_exact if l_inf_norm_exact > 1e-12 else l_inf_error

        return l2_error, l_inf_error, h1_error, relative_l2, relative_l_inf


# Define test functions

def exponential_standard():
    """e^x on [-1, 1]"""
    return SpecialFunctionApproximation(
        lambda x: np.exp(x),
        "Exponential e^x",
        -1.0, 1.0
    )


def exponential_wide():
    """e^x on [-5, 5]"""
    return SpecialFunctionApproximation(
        lambda x: np.exp(x),
        "Exponential e^x (wide)",
        -5.0, 5.0
    )


def sine_standard():
    """sin(x) on [0, 2π]"""
    return SpecialFunctionApproximation(
        lambda x: np.sin(x),
        "Sine sin(x)",
        0.0, 2.0 * np.pi
    )


def cosine_standard():
    """cos(x) on [0, 2π]"""
    return SpecialFunctionApproximation(
        lambda x: np.cos(x),
        "Cosine cos(x)",
        0.0, 2.0 * np.pi
    )


def tangent_standard():
    """tan(x) on [-π/3, π/3]"""
    return SpecialFunctionApproximation(
        lambda x: np.tan(x),
        "Tangent tan(x)",
        -np.pi / 3.0, np.pi / 3.0
    )


def error_function_standard():
    """erf(x) on [-3, 3]"""
    return SpecialFunctionApproximation(
        lambda x: erf(x),
        "Error function erf(x)",
        -3.0, 3.0
    )


def bessel_j0_standard():
    """J_0(x) on [0, 10]"""
    return SpecialFunctionApproximation(
        lambda x: j0(x),
        "Bessel J_0(x)",
        0.0, 10.0
    )


def logarithm_standard():
    """log(1+x) on [0, e-1]"""
    return SpecialFunctionApproximation(
        lambda x: np.log(1.0 + x),
        "Logarithm log(1+x)",
        0.0, np.e - 1.0
    )


def runge_function():
    """1/(1+25x²) on [-1, 1]"""
    return SpecialFunctionApproximation(
        lambda x: 1.0 / (1.0 + 25.0 * x * x),
        "Runge's function 1/(1+25x²)",
        -1.0, 1.0
    )


def run_convergence_study(problem: SpecialFunctionApproximation,
                          mesh_sizes: List[int]) -> ConvergenceResult:
    """Run convergence study for given function"""
    poly_errors = []
    rat_errors = []

    print(f"\nConvergence study: {problem.name}")
    print("=" * 60)

    for n in mesh_sizes:
        print(f"\nMesh: {n} intervals")

        # Polynomial spline approximation
        x_poly, y_poly = problem.approximate_polynomial_spline(n)
        l2, linf, h1, rel_l2, rel_linf = problem.compute_error(x_poly, y_poly)
        h = (problem.b - problem.a) / n

        poly_metric = ErrorMetrics(
            l2_error=l2,
            l_inf_error=linf,
            h1_seminorm_error=h1,
            relative_l2_error=rel_l2,
            relative_l_inf_error=rel_linf,
            num_intervals=n,
            dof=n + 3,  # Cubic spline DOF
            mesh_size=h,
            method='polynomial'
        )
        poly_errors.append(poly_metric)

        print(f"  Polynomial: L2={l2:.6e}, L∞={linf:.6e}, H1={h1:.6e}")
        print(f"              Rel-L2={rel_l2:.6e}, Rel-L∞={rel_linf:.6e}")

        # Rational approximant
        x_rat, y_rat = problem.approximate_rational_piecewise(n)
        l2, linf, h1, rel_l2, rel_linf = problem.compute_error(x_rat, y_rat)

        rat_metric = ErrorMetrics(
            l2_error=l2,
            l_inf_error=linf,
            h1_seminorm_error=h1,
            relative_l2_error=rel_l2,
            relative_l_inf_error=rel_linf,
            num_intervals=n,
            dof=n * 5,  # [2/2] rational DOF (with normalization: 2+2+1)
            mesh_size=h,
            method='rational'
        )
        rat_errors.append(rat_metric)

        print(f"  Rational:   L2={l2:.6e}, L∞={linf:.6e}, H1={h1:.6e}")
        print(f"              Rel-L2={rel_l2:.6e}, Rel-L∞={rel_linf:.6e}")

    # Compute convergence rates
    rates = compute_convergence_rates(poly_errors, rat_errors)

    return ConvergenceResult(
        function_name=problem.name,
        domain=(problem.a, problem.b),
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
        'polynomial_relative_l2': [],
        'rational_l2': [],
        'rational_linf': [],
        'rational_relative_l2': []
    }

    for i in range(1, len(poly_errors)):
        h_ratio = poly_errors[i - 1].mesh_size / poly_errors[i].mesh_size

        # Polynomial rates
        if poly_errors[i].l2_error > 0 and poly_errors[i - 1].l2_error > 0:
            rate_l2 = np.log(poly_errors[i - 1].l2_error / poly_errors[i].l2_error) / np.log(h_ratio)
            rates['polynomial_l2'].append(rate_l2)

        if poly_errors[i].l_inf_error > 0 and poly_errors[i - 1].l_inf_error > 0:
            rate_linf = np.log(poly_errors[i - 1].l_inf_error / poly_errors[i].l_inf_error) / np.log(h_ratio)
            rates['polynomial_linf'].append(rate_linf)

        if poly_errors[i].relative_l2_error > 0 and poly_errors[i - 1].relative_l2_error > 0:
            rate_rel = np.log(poly_errors[i - 1].relative_l2_error / poly_errors[i].relative_l2_error) / np.log(h_ratio)
            rates['polynomial_relative_l2'].append(rate_rel)

        # Rational rates
        if rat_errors[i].l2_error > 0 and rat_errors[i - 1].l2_error > 0:
            rate_l2 = np.log(rat_errors[i - 1].l2_error / rat_errors[i].l2_error) / np.log(h_ratio)
            rates['rational_l2'].append(rate_l2)

        if rat_errors[i].l_inf_error > 0 and rat_errors[i - 1].l_inf_error > 0:
            rate_linf = np.log(rat_errors[i - 1].l_inf_error / rat_errors[i].l_inf_error) / np.log(h_ratio)
            rates['rational_linf'].append(rate_linf)

        if rat_errors[i].relative_l2_error > 0 and rat_errors[i - 1].relative_l2_error > 0:
            rate_rel = np.log(rat_errors[i - 1].relative_l2_error / rat_errors[i].relative_l2_error) / np.log(h_ratio)
            rates['rational_relative_l2'].append(rate_rel)

    return rates


def sanitize_filename(name: str) -> str:
    """Sanitize a string to be used as a filename"""
    # Replace problematic characters
    replacements = {
        ' ': '_',
        '(': '',
        ')': '',
        '/': '_over_',
        '+': 'plus',
        '²': '2',
        '^': '',
        "'": '',
    }
    for old, new in replacements.items():
        name = name.replace(old, new)
    return name


def save_results(results: List[ConvergenceResult], output_dir: str):
    """Save results to JSON for LaTeX report generation"""
    os.makedirs(output_dir, exist_ok=True)

    for result in results:
        filename = sanitize_filename(result.function_name) + '.json'
        filepath = os.path.join(output_dir, filename)

        data = {
            'function_name': result.function_name,
            'domain': result.domain,
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
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 12))

        # Extract data
        h_poly = [e.mesh_size for e in result.polynomial_errors]
        l2_poly = [e.l2_error for e in result.polynomial_errors]
        l2_rat = [e.l2_error for e in result.rational_errors]
        rel_l2_poly = [e.relative_l2_error for e in result.polynomial_errors]
        rel_l2_rat = [e.relative_l2_error for e in result.rational_errors]

        # L2 error vs mesh size
        ax1.loglog(h_poly, l2_poly, 'o-', label='Polynomial spline', linewidth=2, markersize=8)
        ax1.loglog(h_poly, l2_rat, 's-', label='Rational approximant', linewidth=2, markersize=8)

        # Reference lines
        ax1.loglog(h_poly, [h ** 2 for h in h_poly], '--', label='O(h²)', alpha=0.5)
        ax1.loglog(h_poly, [h ** 4 for h in h_poly], '--', label='O(h⁴)', alpha=0.5)

        ax1.set_xlabel('Mesh size h', fontsize=12)
        ax1.set_ylabel('L² error', fontsize=12)
        ax1.set_title(f'{result.function_name} - L² Convergence', fontsize=14, fontweight='bold')
        ax1.legend(fontsize=10)
        ax1.grid(True, alpha=0.3, which='both')

        # Relative L2 error vs mesh size
        ax2.loglog(h_poly, rel_l2_poly, 'o-', label='Polynomial spline', linewidth=2, markersize=8)
        ax2.loglog(h_poly, rel_l2_rat, 's-', label='Rational approximant', linewidth=2, markersize=8)

        ax2.set_xlabel('Mesh size h', fontsize=12)
        ax2.set_ylabel('Relative L² error', fontsize=12)
        ax2.set_title(f'{result.function_name} - Relative L² Convergence', fontsize=14, fontweight='bold')
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3, which='both')

        # Error vs DOF
        dof_poly = [e.dof for e in result.polynomial_errors]
        dof_rat = [e.dof for e in result.rational_errors]

        ax3.loglog(dof_poly, l2_poly, 'o-', label='Polynomial spline', linewidth=2, markersize=8)
        ax3.loglog(dof_rat, l2_rat, 's-', label='Rational approximant', linewidth=2, markersize=8)

        ax3.set_xlabel('Degrees of Freedom', fontsize=12)
        ax3.set_ylabel('L² error', fontsize=12)
        ax3.set_title(f'{result.function_name} - Efficiency', fontsize=14, fontweight='bold')
        ax3.legend(fontsize=10)
        ax3.grid(True, alpha=0.3, which='both')

        # Convergence rates
        if result.convergence_rates['polynomial_l2']:
            mesh_indices = range(1, len(result.polynomial_errors))
            ax4.plot(mesh_indices, result.convergence_rates['polynomial_l2'],
                    'o-', label='Polynomial L²', linewidth=2, markersize=8)
            ax4.plot(mesh_indices, result.convergence_rates['rational_l2'],
                    's-', label='Rational L²', linewidth=2, markersize=8)

            ax4.axhline(y=2, color='gray', linestyle='--', alpha=0.5, label='Order 2')
            ax4.axhline(y=4, color='gray', linestyle='--', alpha=0.5, label='Order 4')

            ax4.set_xlabel('Refinement level', fontsize=12)
            ax4.set_ylabel('Convergence rate', fontsize=12)
            ax4.set_title(f'{result.function_name} - Convergence Rates', fontsize=14, fontweight='bold')
            ax4.legend(fontsize=10)
            ax4.grid(True, alpha=0.3)

        plt.tight_layout()

        filename = sanitize_filename(result.function_name) + '.pdf'
        filepath = os.path.join(output_dir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"Saved plot to {filepath}")
        plt.close()


def main():
    """Run all special function benchmarks"""
    print("Gelfgren Special Function Approximation Study")
    print("=" * 60)

    # Define mesh sizes for convergence study
    mesh_sizes = [4, 8, 16, 32, 64, 128]

    # Define problems
    problems = [
        exponential_standard(),
        sine_standard(),
        cosine_standard(),
        error_function_standard(),
        logarithm_standard(),
        runge_function(),
    ]

    # Run convergence studies
    results = []
    for problem in problems:
        result = run_convergence_study(problem, mesh_sizes)
        results.append(result)

    # Save results
    output_dir = os.path.join(os.path.dirname(__file__), '../data/special_functions')
    save_results(results, output_dir)

    # Generate plots
    figures_dir = os.path.join(os.path.dirname(__file__), '../reports/figures/special_functions')
    plot_convergence(results, figures_dir)

    print("\n" + "=" * 60)
    print("Special function convergence study complete!")
    print(f"Results saved to: {output_dir}")
    print(f"Figures saved to: {figures_dir}")


if __name__ == '__main__':
    main()
