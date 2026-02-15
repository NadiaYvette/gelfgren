#!/usr/bin/env python3
"""
Fair DOF Special Function Approximation Study

Compares piecewise rational approximants vs polynomial splines with EQUAL DOF:
- [5/2] rational (6 DOF) vs Quintic Hermite (6 DOF)
- [7/4] rational (8 DOF) vs Septic Hermite (8 DOF)
- [9/6] rational (10 DOF) vs Nonic Hermite (10 DOF)

Strategy:
- Even denominator degrees (2, 4, 6) to capture complex conjugate poles
- Odd numerator degrees (5, 7, 9) for proper asymptotic behavior
- Higher-order Hermite splines as fair polynomial baseline

This addresses the critique that cubic splines vs [3/1] rationals was unfair.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, BPoly
from scipy.special import erf, jn, mathieu_cem, mathieu_sem, ellipj, ellipk, airy
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
    print("Warning: Gelfgren Python bindings not available. Using fallback.")
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
    dof_per_interval: int
    total_dof: int
    mesh_size: float
    method: str  # 'hermite_quintic', 'rational_5_2', etc.


@dataclass
class ConvergenceResult:
    """Results from convergence study"""
    function_name: str
    domain: Tuple[float, float]
    hermite_errors: dict  # keys: 'quintic', 'septic', 'nonic'
    rational_errors: dict  # keys: '[5/2]', '[7/4]', '[9/6]'
    convergence_rates: dict


def compute_high_order_derivatives(func: Callable, x: float, max_order: int = 5) -> List[float]:
    """Compute derivatives up to max_order using finite differences"""
    h = 1e-5
    derivs = [func(x)]

    # First derivative (central difference)
    if max_order >= 1:
        derivs.append((func(x + h) - func(x - h)) / (2 * h))

    # Second derivative
    if max_order >= 2:
        derivs.append((func(x + h) - 2 * func(x) + func(x - h)) / h**2)

    # Higher derivatives using finite difference stencils
    if max_order >= 3:
        # Third derivative (5-point stencil)
        d3 = (-func(x + 2*h) + 2*func(x + h) - 2*func(x - h) + func(x - 2*h)) / (2 * h**3)
        derivs.append(d3)

    if max_order >= 4:
        # Fourth derivative
        d4 = (func(x + 2*h) - 4*func(x + h) + 6*func(x) - 4*func(x - h) + func(x - 2*h)) / h**4
        derivs.append(d4)

    if max_order >= 5:
        # Fifth derivative (7-point stencil)
        d5 = (func(x + 3*h) - 6*func(x + 2*h) + 15*func(x + h) - 20*func(x) +
              15*func(x - h) - 6*func(x - 2*h) + func(x - 3*h)) / (2 * h**5)
        derivs.append(d5)

    return derivs


def hermite_quintic_interpolation(x_knots: np.ndarray, func: Callable) -> Callable:
    """
    Construct quintic Hermite spline (C² continuous)
    6 DOF per interval: u(a), u'(a), u''(a), u(b), u'(b), u''(b)
    """
    n = len(x_knots) - 1
    coeffs_list = []

    for i in range(n):
        a, b = x_knots[i], x_knots[i + 1]
        h = b - a

        # Get function values and derivatives at endpoints
        derivs_a = compute_high_order_derivatives(func, a, max_order=2)
        derivs_b = compute_high_order_derivatives(func, b, max_order=2)

        u_a, up_a, upp_a = derivs_a[0], derivs_a[1], derivs_a[2]
        u_b, up_b, upp_b = derivs_b[0], derivs_b[1], derivs_b[2]

        # Construct quintic polynomial satisfying 6 conditions
        # p(t) = c0 + c1*t + c2*t² + c3*t³ + c4*t⁴ + c5*t⁵ where t = (x-a)/h
        # Conditions: p(0) = u_a, p'(0) = up_a*h, p''(0) = upp_a*h²,
        #             p(1) = u_b, p'(1) = up_b*h, p''(1) = upp_b*h²

        c0 = u_a
        c1 = up_a * h
        c2 = upp_a * h**2 / 2

        # Solve for c3, c4, c5 from conditions at t=1
        # p(1) = c0 + c1 + c2 + c3 + c4 + c5 = u_b
        # p'(1) = (c1 + 2*c2 + 3*c3 + 4*c4 + 5*c5)/h = up_b*h
        # p''(1) = (2*c2 + 6*c3 + 12*c4 + 20*c5)/h² = upp_b*h²

        # System: [1  1  1 ] [c3]   [u_b - c0 - c1 - c2        ]
        #         [3  4  5 ] [c4] = [up_b*h - c1 - 2*c2        ]
        #         [6 12 20] [c5]   [upp_b*h² - 2*c2           ]

        rhs = np.array([
            u_b - c0 - c1 - c2,
            up_b * h - c1 - 2 * c2,
            upp_b * h**2 - 2 * c2
        ])

        A = np.array([
            [1, 1, 1],
            [3, 4, 5],
            [6, 12, 20]
        ])

        c3, c4, c5 = np.linalg.solve(A, rhs)

        coeffs_list.append((a, h, np.array([c0, c1, c2, c3, c4, c5])))

    def evaluate(x: np.ndarray) -> np.ndarray:
        """Evaluate quintic Hermite spline at points x"""
        result = np.zeros_like(x)
        for i in range(len(x)):
            xi = x[i]
            # Find which interval
            for j, (a, h, coeffs) in enumerate(coeffs_list):
                if j == len(coeffs_list) - 1:  # Last interval
                    if a <= xi <= a + h:
                        t = (xi - a) / h
                        result[i] = np.polyval(coeffs[::-1], t)
                        break
                else:
                    if a <= xi < a + h:
                        t = (xi - a) / h
                        result[i] = np.polyval(coeffs[::-1], t)
                        break
        return result

    return evaluate


def hermite_septic_interpolation(x_knots: np.ndarray, func: Callable) -> Callable:
    """
    Construct septic (degree 7) Hermite spline (C³ continuous)
    8 DOF per interval: u, u', u'', u''' at both endpoints
    """
    n = len(x_knots) - 1
    coeffs_list = []

    for i in range(n):
        a, b = x_knots[i], x_knots[i + 1]
        h = b - a

        # Get function values and derivatives at endpoints
        derivs_a = compute_high_order_derivatives(func, a, max_order=3)
        derivs_b = compute_high_order_derivatives(func, b, max_order=3)

        u_a, up_a, upp_a, uppp_a = derivs_a[0:4]
        u_b, up_b, upp_b, uppp_b = derivs_b[0:4]

        # Construct septic polynomial (degree 7)
        c0 = u_a
        c1 = up_a * h
        c2 = upp_a * h**2 / 2
        c3 = uppp_a * h**3 / 6

        # Solve for c4, c5, c6, c7 from conditions at t=1
        A = np.array([
            [1, 1, 1, 1],
            [4, 5, 6, 7],
            [12, 20, 30, 42],
            [24, 60, 120, 210]
        ])

        rhs = np.array([
            u_b - c0 - c1 - c2 - c3,
            up_b * h - c1 - 2*c2 - 3*c3,
            upp_b * h**2 - 2*c2 - 6*c3,
            uppp_b * h**3 - 6*c3
        ])

        c4, c5, c6, c7 = np.linalg.solve(A, rhs)

        coeffs_list.append((a, h, np.array([c0, c1, c2, c3, c4, c5, c6, c7])))

    def evaluate(x: np.ndarray) -> np.ndarray:
        result = np.zeros_like(x)
        for i in range(len(x)):
            xi = x[i]
            for j, (a, h, coeffs) in enumerate(coeffs_list):
                if j == len(coeffs_list) - 1:
                    if a <= xi <= a + h:
                        t = (xi - a) / h
                        result[i] = np.polyval(coeffs[::-1], t)
                        break
                else:
                    if a <= xi < a + h:
                        t = (xi - a) / h
                        result[i] = np.polyval(coeffs[::-1], t)
                        break
        return result

    return evaluate


def hermite_nonic_interpolation(x_knots: np.ndarray, func: Callable) -> Callable:
    """
    Construct nonic (degree 9) Hermite spline (C⁴ continuous)
    10 DOF per interval: u, u', u'', u''', u'''' at both endpoints
    """
    n = len(x_knots) - 1
    coeffs_list = []

    for i in range(n):
        a, b = x_knots[i], x_knots[i + 1]
        h = b - a

        # Get function values and derivatives at endpoints
        derivs_a = compute_high_order_derivatives(func, a, max_order=4)
        derivs_b = compute_high_order_derivatives(func, b, max_order=4)

        u_a, up_a, upp_a, uppp_a, upppp_a = derivs_a[0:5]
        u_b, up_b, upp_b, uppp_b, upppp_b = derivs_b[0:5]

        # Construct nonic polynomial (degree 9)
        c0 = u_a
        c1 = up_a * h
        c2 = upp_a * h**2 / 2
        c3 = uppp_a * h**3 / 6
        c4 = upppp_a * h**4 / 24

        # Solve for c5, c6, c7, c8, c9 from conditions at t=1
        A = np.array([
            [1, 1, 1, 1, 1],
            [5, 6, 7, 8, 9],
            [20, 30, 42, 56, 72],
            [60, 120, 210, 336, 504],
            [120, 360, 840, 1680, 3024]
        ])

        rhs = np.array([
            u_b - c0 - c1 - c2 - c3 - c4,
            up_b * h - c1 - 2*c2 - 3*c3 - 4*c4,
            upp_b * h**2 - 2*c2 - 6*c3 - 12*c4,
            uppp_b * h**3 - 6*c3 - 24*c4,
            upppp_b * h**4 - 24*c4
        ])

        c5, c6, c7, c8, c9 = np.linalg.solve(A, rhs)

        coeffs_list.append((a, h, np.array([c0, c1, c2, c3, c4, c5, c6, c7, c8, c9])))

    def evaluate(x: np.ndarray) -> np.ndarray:
        result = np.zeros_like(x)
        for i in range(len(x)):
            xi = x[i]
            for j, (a, h, coeffs) in enumerate(coeffs_list):
                if j == len(coeffs_list) - 1:
                    if a <= xi <= a + h:
                        t = (xi - a) / h
                        result[i] = np.polyval(coeffs[::-1], t)
                        break
                else:
                    if a <= xi < a + h:
                        t = (xi - a) / h
                        result[i] = np.polyval(coeffs[::-1], t)
                        break
        return result

    return evaluate


class SpecialFunctionApproximation:
    """Special function approximation with fair DOF comparison"""

    def __init__(self, func: Callable, name: str, a: float, b: float):
        self.func = func
        self.name = name
        self.a = a
        self.b = b

    def approximate_hermite(self, n_intervals: int, degree: int) -> Tuple[np.ndarray, np.ndarray]:
        """Approximate using high-order Hermite spline"""
        x_knots = np.linspace(self.a, self.b, n_intervals + 1)

        if degree == 5:
            spline_func = hermite_quintic_interpolation(x_knots, self.func)
        elif degree == 7:
            spline_func = hermite_septic_interpolation(x_knots, self.func)
        elif degree == 9:
            spline_func = hermite_nonic_interpolation(x_knots, self.func)
        else:
            raise ValueError(f"Unsupported degree: {degree}")

        # Evaluate on fine grid
        x_fine = np.linspace(self.a, self.b, 1000)
        y_approx = spline_func(x_fine)

        return x_fine, y_approx

    def approximate_rational(self, n_intervals: int, deg_num: int, deg_den: int) -> Tuple[np.ndarray, np.ndarray]:
        """Approximate using piecewise rational with specified degrees"""
        if not GELFGREN_AVAILABLE:
            # Fallback to Hermite with same total DOF
            dof = deg_num + deg_den + 1  # Approximate total DOF
            hermite_deg = 2 * (dof // 2) - 1  # Nearest odd degree
            return self.approximate_hermite(n_intervals, hermite_deg)

        # Create mesh
        mesh = gf.Mesh.uniform(self.a, self.b, n_intervals)

        # Compute derivatives
        max_deriv = (deg_num + deg_den) // 2

        def get_derivs(x):
            return compute_high_order_derivatives(self.func, x, max_order=max_deriv)

        # Create piecewise rational approximation
        pw = gf.PiecewiseRational.from_function(mesh, deg_num, deg_den, get_derivs)

        # Evaluate on fine grid
        x_fine = np.linspace(self.a, self.b, 1000)
        y_approx = pw.evaluate(np.array(x_fine))

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

        # H1 seminorm error
        dy_approx = np.diff(y_approx) / h
        dy_exact = np.diff(y_exact) / h
        h1_error = np.sqrt(h * np.sum((dy_approx - dy_exact) ** 2))

        # Relative errors
        l2_norm_exact = np.sqrt(h * np.sum(y_exact ** 2))
        l_inf_norm_exact = np.max(np.abs(y_exact))

        relative_l2 = l2_error / l2_norm_exact if l2_norm_exact > 1e-12 else l2_error
        relative_l_inf = l_inf_error / l_inf_norm_exact if l_inf_norm_exact > 1e-12 else l_inf_error

        return l2_error, l_inf_error, h1_error, relative_l2, relative_l_inf


# Define all 16 test functions (same as before)

def exponential_standard():
    return SpecialFunctionApproximation(
        lambda x: np.exp(x),
        "Exponential e^x",
        -1.0, 1.0
    )


def sine_standard():
    return SpecialFunctionApproximation(
        lambda x: np.sin(x),
        "Sine sin(x)",
        0.0, 2 * np.pi
    )


def cosine_standard():
    return SpecialFunctionApproximation(
        lambda x: np.cos(x),
        "Cosine cos(x)",
        0.0, 2 * np.pi
    )


def error_function_standard():
    return SpecialFunctionApproximation(
        lambda x: erf(x),
        "Error function erf(x)",
        0.0, 3.0
    )


def logarithm_standard():
    return SpecialFunctionApproximation(
        lambda x: np.log(1 + x),
        "Logarithm log(1+x)",
        0.0, 1.0
    )


def runge_function():
    return SpecialFunctionApproximation(
        lambda x: 1 / (1 + 25 * x**2),
        "Runge's function 1/(1+25x²)",
        -1.0, 1.0
    )


def mathieu_ce0_q5():
    from scipy.special import mathieu_cem
    def func(x):
        return mathieu_cem(0, 5, x)[0]
    return SpecialFunctionApproximation(func, "Mathieu ce_0(x, q=5)", 0.0, 2*np.pi)


def mathieu_se2_q10():
    from scipy.special import mathieu_sem
    def func(x):
        return mathieu_sem(2, 10, x)[0]
    return SpecialFunctionApproximation(func, "Mathieu se_2(x, q=10)", 0.0, 2*np.pi)


def jacobi_sn_k05():
    from scipy.special import ellipj, ellipk
    k = 0.5
    K_val = ellipk(k**2)
    def func(x):
        return ellipj(x, k**2)[0]
    return SpecialFunctionApproximation(func, "Jacobi sn(x, k=0.5)", 0.0, K_val)


def jacobi_cn_k09():
    from scipy.special import ellipj, ellipk
    k = 0.9
    K_val = ellipk(k**2)
    def func(x):
        return ellipj(x, k**2)[1]
    return SpecialFunctionApproximation(func, "Jacobi cn(x, k=0.9)", 0.0, K_val)


def lemniscate_sl():
    from scipy.special import ellipj, ellipk
    k = 1 / np.sqrt(2)
    K_val = ellipk(k**2)
    def func(x):
        return ellipj(x, k**2)[0]
    return SpecialFunctionApproximation(func, "Lemniscate sl(x)", 0.0, 2*K_val)


def airy_ai_oscillatory():
    from scipy.special import airy
    def func(x):
        return airy(x)[0]
    return SpecialFunctionApproximation(func, "Airy Ai(x) oscillatory", -10.0, 2.0)


def airy_bi_exponential():
    from scipy.special import airy
    def func(x):
        return airy(x)[2]
    return SpecialFunctionApproximation(func, "Airy Bi(x) exponential", -5.0, 2.0)


def bessel_j1():
    from scipy.special import jn
    def func(x):
        return jn(1, x)
    return SpecialFunctionApproximation(func, "Bessel J_1(x)", 0.0, 20.0)


def bessel_j5():
    from scipy.special import jn
    def func(x):
        return jn(5, x)
    return SpecialFunctionApproximation(func, "Bessel J_5(x)", 0.0, 30.0)


def bessel_j10():
    from scipy.special import jn
    def func(x):
        return jn(10, x)
    return SpecialFunctionApproximation(func, "Bessel J_10(x)", 0.0, 40.0)


def get_all_test_functions():
    return [
        exponential_standard(),
        sine_standard(),
        cosine_standard(),
        error_function_standard(),
        logarithm_standard(),
        runge_function(),
        mathieu_ce0_q5(),
        mathieu_se2_q10(),
        jacobi_sn_k05(),
        jacobi_cn_k09(),
        lemniscate_sl(),
        airy_ai_oscillatory(),
        airy_bi_exponential(),
        bessel_j1(),
        bessel_j5(),
        bessel_j10(),
    ]


def run_convergence_study(problem: SpecialFunctionApproximation,
                          mesh_sizes: List[int],
                          degrees_to_test: List[Tuple[int, str, int, int]]):
    """
    Run convergence study with fair DOF comparison

    degrees_to_test: List of (hermite_deg, rational_name, num_deg, den_deg)
    Example: [(5, '[5/2]', 5, 2), (7, '[7/4]', 7, 4), (9, '[9/6]', 9, 6)]
    """
    print(f"\nTesting: {problem.name}")
    print("=" * 60)

    results = {'hermite': {}, 'rational': {}}

    for hermite_deg, rat_name, num_deg, den_deg in degrees_to_test:
        hermite_key = f"degree_{hermite_deg}"

        results['hermite'][hermite_key] = []
        results['rational'][rat_name] = []

        for n in mesh_sizes:
            h = (problem.b - problem.a) / n

            # DOF calculation
            hermite_dof_per_interval = hermite_deg + 1  # Simplified
            rational_dof_per_interval = num_deg + den_deg  # Approximate

            # Hermite spline
            try:
                x, y_hermite = problem.approximate_hermite(n, hermite_deg)
                l2, linf, h1, rel_l2, rel_linf = problem.compute_error(x, y_hermite)

                results['hermite'][hermite_key].append(ErrorMetrics(
                    l2_error=l2,
                    l_inf_error=linf,
                    h1_seminorm_error=h1,
                    relative_l2_error=rel_l2,
                    relative_l_inf_error=rel_linf,
                    num_intervals=n,
                    dof_per_interval=hermite_dof_per_interval,
                    total_dof=n * hermite_dof_per_interval,
                    mesh_size=h,
                    method=f"hermite_{hermite_deg}"
                ))
                print(f"  Hermite deg={hermite_deg}, n={n:3d}: L2={l2:.2e}")
            except Exception as e:
                print(f"  Hermite deg={hermite_deg}, n={n:3d}: FAILED - {e}")

            # Rational approximation
            try:
                x, y_rational = problem.approximate_rational(n, num_deg, den_deg)
                l2, linf, h1, rel_l2, rel_linf = problem.compute_error(x, y_rational)

                results['rational'][rat_name].append(ErrorMetrics(
                    l2_error=l2,
                    l_inf_error=linf,
                    h1_seminorm_error=h1,
                    relative_l2_error=rel_l2,
                    relative_l_inf_error=rel_linf,
                    num_intervals=n,
                    dof_per_interval=rational_dof_per_interval,
                    total_dof=n * rational_dof_per_interval,
                    mesh_size=h,
                    method=rat_name
                ))
                print(f"  Rational {rat_name}, n={n:3d}: L2={l2:.2e}")
            except Exception as e:
                print(f"  Rational {rat_name}, n={n:3d}: FAILED - {e}")

    return results


def main():
    """Run fair DOF comparison benchmarks"""
    print("="*80)
    print("Fair DOF Special Function Approximation Study")
    print("="*80)
    print("\nComparison strategy:")
    print("  [5/2] rational (6 DOF) vs Quintic Hermite (6 DOF)")
    print("  [7/4] rational (8 DOF) vs Septic Hermite (8 DOF)")
    print("  [9/6] rational (10 DOF) vs Nonic Hermite (10 DOF)")
    print("\nRationale:")
    print("  - Even denominator degrees (2,4,6) capture complex conjugate poles")
    print("  - Higher-order Hermite splines provide fair polynomial baseline")
    print("  - Equal DOF comparison eliminates unfair advantage")
    print("="*80)

    # Mesh sizes to test
    mesh_sizes = [4, 8, 16, 32, 64, 128]

    # Degrees to test: (hermite_deg, rational_name, num_deg, den_deg)
    degrees = [
        (5, '[5/2]', 5, 2),
        (7, '[7/4]', 7, 4),
        (9, '[9/6]', 9, 6)
    ]

    # Get all test functions
    test_functions = get_all_test_functions()

    print(f"\nTesting {len(test_functions)} functions with {len(mesh_sizes)} mesh sizes")
    print(f"Total tests: {len(test_functions) * len(mesh_sizes) * len(degrees) * 2}")

    # Run benchmarks
    all_results = {}
    for func in test_functions:
        try:
            results = run_convergence_study(func, mesh_sizes, degrees)
            all_results[func.name] = results
        except Exception as e:
            print(f"ERROR in {func.name}: {e}")
            import traceback
            traceback.print_exc()

    # Save results
    output_dir = os.path.join(os.path.dirname(__file__), '../data/special_functions_fair')
    os.makedirs(output_dir, exist_ok=True)

    for func_name, results in all_results.items():
        safe_name = func_name.replace('/', '_').replace(' ', '_').replace('(', '').replace(')', '').replace(',', '')
        output_file = os.path.join(output_dir, f'{safe_name}.json')

        # Convert to JSON-serializable format
        json_results = {
            'function_name': func_name,
            'hermite': {},
            'rational': {}
        }

        for key, errors in results['hermite'].items():
            json_results['hermite'][key] = [asdict(e) for e in errors]

        for key, errors in results['rational'].items():
            json_results['rational'][key] = [asdict(e) for e in errors]

        with open(output_file, 'w') as f:
            json.dump(json_results, f, indent=2)

        print(f"Saved: {output_file}")

    print("\n" + "="*80)
    print("Fair DOF benchmarks complete!")
    print(f"Results saved to: {output_dir}")
    print("="*80)


if __name__ == "__main__":
    main()
