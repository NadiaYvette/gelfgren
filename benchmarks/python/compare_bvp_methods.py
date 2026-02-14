#!/usr/bin/env python3
"""
Comprehensive BVP Method Comparison

Compares three approaches:
1. Polynomial Finite Differences (baseline)
2. Rational Collocation - Quadratic Form
3. Rational Collocation - Cleared Form
"""

import numpy as np
import json
from pathlib import Path
from dataclasses import dataclass, asdict
import time

from rational_collocation import RationalCollocationQuadratic, QConstraintType
from rational_collocation_cleared import RationalCollocationCleared


@dataclass
class MethodResult:
    """Results from a single method on a single problem"""
    problem_name: str
    method_name: str
    formulation: str  # "Polynomial FD", "Quadratic", "Cleared"
    degree_n: int
    degree_m: int
    n_collocation: int
    max_error: float
    l2_error: float
    solve_time_ms: float
    success: bool
    residual_norm: float = 0.0


class PolynomialFiniteDifference:
    """Polynomial finite difference BVP solver (baseline)"""

    def __init__(self, f, a, b, alpha, beta):
        self.f = f
        self.a = a
        self.b = b
        self.alpha = alpha
        self.beta = beta

    def solve(self, n):
        """Solve with n interior points"""
        x = np.linspace(self.a, self.b, n + 2)
        h = x[1] - x[0]

        A = np.zeros((n + 2, n + 2))
        b = np.zeros(n + 2)

        for i in range(1, n + 1):
            A[i, i-1] = -1.0 / h**2
            A[i, i] = 2.0 / h**2
            A[i, i+1] = -1.0 / h**2
            b[i] = self.f(x[i])

        A[0, 0] = 1.0
        A[n + 1, n + 1] = 1.0
        b[0] = self.alpha
        b[n + 1] = self.beta

        u = np.linalg.solve(A, b)
        return x, u, True


def test_problem_smooth_poisson():
    """u(x) = x(1-x), -u'' = 2"""
    return {
        'name': 'Smooth_Poisson',
        'f': lambda x: 2.0,
        'exact': lambda x: x * (1 - x),
        'a': 0.0,
        'b': 1.0,
        'alpha': 0.0,
        'beta': 0.0,
    }


def test_problem_smooth_trig():
    """u(x) = sin(πx), -u'' = π² sin(πx)"""
    return {
        'name': 'Smooth_Trig',
        'f': lambda x: np.pi**2 * np.sin(np.pi * x),
        'exact': lambda x: np.sin(np.pi * x),
        'a': 0.0,
        'b': 1.0,
        'alpha': 0.0,
        'beta': 0.0,
    }


def test_problem_smooth_exp():
    """Exponential solution"""
    alpha = 1.0
    c1 = -alpha / (np.exp(1) - 1)
    c2 = alpha

    def exact(x):
        return c1 * np.exp(x) + c2 * np.exp(-x)

    def f(x):
        return -c1 * np.exp(x) - c2 * np.exp(-x)

    return {
        'name': 'Smooth_Exp',
        'f': f,
        'exact': exact,
        'a': 0.0,
        'b': 1.0,
        'alpha': 0.0,
        'beta': 0.0,
    }


def run_comparison(problems, poly_grid_sizes, rational_degrees):
    """
    Run comprehensive comparison of all three methods.

    Parameters
    ----------
    problems : list of dict
        Test problems
    poly_grid_sizes : list of int
        Grid sizes for polynomial FD (n interior points)
    rational_degrees : list of (n, m) tuples
        Rational approximant degrees to test
    """
    all_results = []

    for problem in problems:
        print(f"\n{'='*70}")
        print(f"Problem: {problem['name']}")
        print(f"{'='*70}\n")

        f = problem['f']
        exact = problem['exact']
        a, b = problem['a'], problem['b']
        alpha, beta = problem['alpha'], problem['beta']

        # Evaluation points for error computation
        x_eval = np.linspace(a, b, 1000)
        u_exact = exact(x_eval)

        # Test polynomial finite differences
        print(f"{'Method':<30} {'DOF':<8} {'Max Error':<15} {'L2 Error':<15} {'Time (ms)':<12}")
        print(f"{'-'*80}")

        poly_solver = PolynomialFiniteDifference(f, a, b, alpha, beta)

        for n_grid in poly_grid_sizes:
            start_time = time.time()
            x_poly, u_poly, success = poly_solver.solve(n_grid)
            solve_time = time.time() - start_time

            u_computed = np.interp(x_eval, x_poly, u_poly)
            error = np.abs(u_computed - u_exact)
            max_error = np.max(error)
            l2_error = np.sqrt(np.mean(error**2))

            print(f"{'Polynomial FD':<30} {n_grid:<8} {max_error:<15.2e} {l2_error:<15.2e} {solve_time*1000:<12.2f}")

            all_results.append(MethodResult(
                problem_name=problem['name'],
                method_name='Polynomial_FD',
                formulation='Polynomial FD',
                degree_n=n_grid,
                degree_m=0,
                n_collocation=n_grid,
                max_error=max_error,
                l2_error=l2_error,
                solve_time_ms=solve_time * 1000,
                success=success
            ))

        # Test rational collocation methods
        for n, m in rational_degrees:
            # Quadratic form
            try:
                solver_quad = RationalCollocationQuadratic(
                    f=f, a=a, b=b, alpha=alpha, beta=beta,
                    n=n, m=m,
                    q_constraint=QConstraintType.ENDPOINT,
                    constraint_epsilon=1e-4
                )

                start_time = time.time()
                result_quad = solver_quad.solve(method='trf', verbose=False)
                solve_time_quad = time.time() - start_time

                u_computed_quad = solver_quad.evaluate_solution(result_quad, x_eval)
                error_quad = np.abs(u_computed_quad - u_exact)
                max_error_quad = np.max(error_quad)
                l2_error_quad = np.sqrt(np.mean(error_quad**2))

                print(f"{'Rational Quadratic [' + str(n) + '/' + str(m) + ']':<30} {solver_quad.k:<8} "
                      f"{max_error_quad:<15.2e} {l2_error_quad:<15.2e} {solve_time_quad*1000:<12.2f}")

                all_results.append(MethodResult(
                    problem_name=problem['name'],
                    method_name='Rational_Quadratic',
                    formulation='Quadratic',
                    degree_n=n,
                    degree_m=m,
                    n_collocation=solver_quad.k,
                    max_error=max_error_quad,
                    l2_error=l2_error_quad,
                    solve_time_ms=solve_time_quad * 1000,
                    success=result_quad.success,
                    residual_norm=result_quad.residual_norm
                ))

            except Exception as e:
                print(f"{'Rational Quadratic [' + str(n) + '/' + str(m) + ']':<30} {'FAILED':<8} {str(e)[:40]}")

            # Cleared form
            try:
                solver_cleared = RationalCollocationCleared(
                    f=f, a=a, b=b, alpha=alpha, beta=beta,
                    n=n, m=m,
                    q_constraint=QConstraintType.ENDPOINT,
                    constraint_epsilon=1e-4
                )

                start_time = time.time()
                result_cleared = solver_cleared.solve(method='trf', verbose=False)
                solve_time_cleared = time.time() - start_time

                u_computed_cleared = solver_cleared.evaluate_solution(result_cleared, x_eval)
                error_cleared = np.abs(u_computed_cleared - u_exact)
                max_error_cleared = np.max(error_cleared)
                l2_error_cleared = np.sqrt(np.mean(error_cleared**2))

                print(f"{'Rational Cleared [' + str(n) + '/' + str(m) + ']':<30} {solver_cleared.k:<8} "
                      f"{max_error_cleared:<15.2e} {l2_error_cleared:<15.2e} {solve_time_cleared*1000:<12.2f}")

                all_results.append(MethodResult(
                    problem_name=problem['name'],
                    method_name='Rational_Cleared',
                    formulation='Cleared',
                    degree_n=n,
                    degree_m=m,
                    n_collocation=solver_cleared.k,
                    max_error=max_error_cleared,
                    l2_error=l2_error_cleared,
                    solve_time_ms=solve_time_cleared * 1000,
                    success=result_cleared.success,
                    residual_norm=result_cleared.residual_norm
                ))

            except Exception as e:
                print(f"{'Rational Cleared [' + str(n) + '/' + str(m) + ']':<30} {'FAILED':<8} {str(e)[:40]}")

    return all_results


def main():
    """Run comprehensive comparison and save results"""

    problems = [
        test_problem_smooth_poisson(),
        test_problem_smooth_trig(),
        test_problem_smooth_exp(),
    ]

    poly_grid_sizes = [10, 20, 40, 80, 160]
    rational_degrees = [(4, 2), (6, 3), (8, 4)]

    print("="*70)
    print("BVP Method Comparison: Polynomial vs Quadratic vs Cleared Forms")
    print("="*70)

    results = run_comparison(problems, poly_grid_sizes, rational_degrees)

    # Save to JSON
    output_dir = Path(__file__).parent.parent / 'data'
    output_file = output_dir / 'bvp_method_comparison.json'

    results_dict = {
        'metadata': {
            'description': 'Comparison of polynomial FD vs rational collocation (quadratic and cleared forms)',
            'polynomial_grid_sizes': poly_grid_sizes,
            'rational_degrees': [f"[{n}/{m}]" for n, m in rational_degrees],
            'constraint_type': 'ENDPOINT',
        },
        'results': [asdict(r) for r in results]
    }

    with open(output_file, 'w') as f:
        json.dump(results_dict, f, indent=2)

    print(f"\n{'='*70}")
    print(f"Results saved to: {output_file}")
    print(f"Total tests: {len(results)}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
