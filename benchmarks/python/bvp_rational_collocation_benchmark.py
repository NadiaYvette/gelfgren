#!/usr/bin/env python3
"""
Benchmark rational collocation BVP solver against polynomial finite differences.

Compares the new quadratic formulation rational collocation with traditional
polynomial finite difference methods on various test problems.
"""

import numpy as np
import json
from pathlib import Path
from rational_collocation import RationalCollocationQuadratic
from dataclasses import dataclass, asdict
from typing import Callable, List, Tuple
import time


@dataclass
class BenchmarkResult:
    """Results from a single benchmark run"""
    problem_name: str
    method: str
    degree: int  # For polynomial or rational degree [n/m]
    n_collocation: int  # Number of collocation points (for rational)
    n_grid: int  # Grid size (for finite difference)
    max_error: float
    l2_error: float
    linf_error: float
    solve_time: float
    success: bool


class PolynomialFiniteDifference:
    """
    Traditional polynomial finite difference BVP solver for comparison.

    Solves: -u'' = f(x) on [a,b] with u(a) = alpha, u(b) = beta
    """

    def __init__(self, f: Callable[[float], float],
                 a: float, b: float,
                 alpha: float, beta: float):
        self.f = f
        self.a = a
        self.b = b
        self.alpha = alpha
        self.beta = beta

    def solve(self, n: int) -> Tuple[np.ndarray, np.ndarray, bool]:
        """
        Solve using finite differences with n interior points.

        Returns
        -------
        x : ndarray
            Grid points
        u : ndarray
            Solution values
        success : bool
            Whether solve succeeded
        """
        # Grid with n interior points
        x = np.linspace(self.a, self.b, n + 2)
        h = x[1] - x[0]

        # Tridiagonal system for -u'' = f
        A = np.zeros((n + 2, n + 2))
        b = np.zeros(n + 2)

        # Interior points: -u''(x_i) = f(x_i)
        # Discretized as: -(u_{i-1} - 2u_i + u_{i+1})/h^2 = f(x_i)
        for i in range(1, n + 1):
            A[i, i-1] = -1.0 / h**2
            A[i, i] = 2.0 / h**2
            A[i, i+1] = -1.0 / h**2
            b[i] = self.f(x[i])

        # Boundary conditions
        A[0, 0] = 1.0
        A[n + 1, n + 1] = 1.0
        b[0] = self.alpha
        b[n + 1] = self.beta

        try:
            u = np.linalg.solve(A, b)
            return x, u, True
        except np.linalg.LinAlgError:
            return x, np.zeros_like(x), False


class BVPBenchmarkSuite:
    """Suite of test problems for BVP benchmarking"""

    @staticmethod
    def smooth_poisson():
        """
        Smooth Poisson problem with polynomial solution.
        -u'' = 2, u(0) = 0, u(1) = 0
        Exact: u(x) = x(1-x)
        """
        return {
            'name': 'Smooth_Poisson',
            'f': lambda x: 2.0,
            'exact': lambda x: x * (1 - x),
            'a': 0.0,
            'b': 1.0,
            'alpha': 0.0,
            'beta': 0.0,
            'description': 'Polynomial solution: u(x) = x(1-x)'
        }

    @staticmethod
    def smooth_polynomial_quartic():
        """
        Smooth problem with quartic polynomial solution.
        -u'' = 12x^2 - 2, u(0) = 0, u(1) = 0
        Exact: u(x) = x^4 - x^2
        """
        return {
            'name': 'Smooth_Quartic',
            'f': lambda x: 12 * x**2 - 2,
            'exact': lambda x: x**4 - x**2,
            'a': 0.0,
            'b': 1.0,
            'alpha': 0.0,
            'beta': 0.0,
            'description': 'Quartic polynomial solution'
        }

    @staticmethod
    def smooth_trig():
        """
        Smooth trigonometric solution.
        -u'' = π^2 sin(πx), u(0) = 0, u(1) = 0
        Exact: u(x) = sin(πx)
        """
        return {
            'name': 'Smooth_Trig',
            'f': lambda x: np.pi**2 * np.sin(np.pi * x),
            'exact': lambda x: np.sin(np.pi * x),
            'a': 0.0,
            'b': 1.0,
            'alpha': 0.0,
            'beta': 0.0,
            'description': 'Trigonometric solution: u(x) = sin(πx)'
        }

    @staticmethod
    def smooth_exp():
        """
        Smooth exponential solution.
        """
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
            'description': 'Exponential solution'
        }


def benchmark_problem(problem: dict, degrees: List[int], grid_sizes: List[int]) -> List[BenchmarkResult]:
    """
    Benchmark a single problem with both methods.

    Parameters
    ----------
    problem : dict
        Problem specification
    degrees : list of int
        Rational degrees to test (will test [n/(n//2)] rationals)
    grid_sizes : list of int
        Grid sizes for finite difference (n interior points)

    Returns
    -------
    list of BenchmarkResult
        Results for all configurations
    """
    results = []

    # Extract problem data
    f = problem['f']
    exact = problem['exact']
    a, b = problem['a'], problem['b']
    alpha, beta = problem['alpha'], problem['beta']
    name = problem['name']

    # Evaluation points for error computation
    x_eval = np.linspace(a, b, 1000)
    u_exact = exact(x_eval)

    print(f"\n{'='*70}")
    print(f"Problem: {name}")
    print(f"{problem['description']}")
    print(f"{'='*70}")

    # Test polynomial finite differences
    print(f"\n{'Polynomial Finite Difference':^70}")
    print(f"{'-'*70}")
    print(f"{'n_grid':<15} {'Max Error':<15} {'L2 Error':<15} {'Time (ms)':<15}")
    print(f"{'-'*70}")

    poly_solver = PolynomialFiniteDifference(f, a, b, alpha, beta)

    for n_grid in grid_sizes:
        start_time = time.time()
        x_poly, u_poly, success = poly_solver.solve(n_grid)
        solve_time = time.time() - start_time

        # Interpolate to evaluation points
        u_computed = np.interp(x_eval, x_poly, u_poly)
        error = np.abs(u_computed - u_exact)

        max_error = np.max(error)
        l2_error = np.sqrt(np.mean(error**2))
        linf_error = max_error

        print(f"{n_grid:<15} {max_error:<15.2e} {l2_error:<15.2e} {solve_time*1000:<15.2f}")

        results.append(BenchmarkResult(
            problem_name=name,
            method='Polynomial_FD',
            degree=n_grid,
            n_collocation=0,
            n_grid=n_grid,
            max_error=max_error,
            l2_error=l2_error,
            linf_error=linf_error,
            solve_time=solve_time,
            success=success
        ))

    # Test rational collocation
    print(f"\n{'Rational Collocation (Quadratic Form)':^70}")
    print(f"{'-'*70}")
    print(f"{'[n/m]':<10} {'k_pts':<10} {'Max Error':<15} {'L2 Error':<15} {'Time (ms)':<15}")
    print(f"{'-'*70}")

    for n in degrees:
        # Use m = n // 2 (e.g., [4/2], [6/3], [8/4])
        m = max(1, n // 2)

        try:
            # Use endpoint constraints to prevent boundary poles
            # while allowing Q to vary in interior (true rational behavior)
            from rational_collocation import QConstraintType

            solver = RationalCollocationQuadratic(
                f=f, a=a, b=b, alpha=alpha, beta=beta,
                n=n, m=m,
                q_constraint=QConstraintType.ENDPOINT,
                constraint_epsilon=1e-4
            )

            start_time = time.time()
            result = solver.solve(method='lm', verbose=False)
            solve_time = time.time() - start_time

            # Check for spurious poles
            has_pole, min_q = check_for_poles(solver, result, a, b)

            if has_pole:
                # Mark as failure - spurious pole detected
                print(f"[{n}/{m}]{'':<5} {solver.k:<10} {'POLE':<15} {f'(Q_min={min_q:.2e})':<15} {solve_time*1000:<15.2f}")
                results.append(BenchmarkResult(
                    problem_name=name,
                    method='Rational_Collocation',
                    degree=n,
                    n_collocation=solver.k,
                    n_grid=0,
                    max_error=float('inf'),
                    l2_error=float('inf'),
                    linf_error=float('inf'),
                    solve_time=solve_time,
                    success=False
                ))
            else:
                # Evaluate at test points
                u_computed = solver.evaluate_solution(result, x_eval)
                error = np.abs(u_computed - u_exact)

                max_error = np.max(error)
                l2_error = np.sqrt(np.mean(error**2))
                linf_error = max_error

                print(f"[{n}/{m}]{'':<5} {solver.k:<10} {max_error:<15.2e} {l2_error:<15.2e} {solve_time*1000:<15.2f}")

                results.append(BenchmarkResult(
                    problem_name=name,
                    method='Rational_Collocation',
                    degree=n,
                    n_collocation=solver.k,
                    n_grid=0,
                    max_error=max_error,
                    l2_error=l2_error,
                    linf_error=linf_error,
                    solve_time=solve_time,
                    success=result.success
                ))

        except Exception as e:
            print(f"[{n}/{m}]{'':<5} {'FAILED':<10} {str(e)[:40]:<40}")
            results.append(BenchmarkResult(
                problem_name=name,
                method='Rational_Collocation',
                degree=n,
                n_collocation=0,
                n_grid=0,
                max_error=float('inf'),
                l2_error=float('inf'),
                linf_error=float('inf'),
                solve_time=0.0,
                success=False
            ))

    return results


def check_for_poles(solver, result, a, b, relative_threshold=1e-6):
    """
    Check if Q has near-zero values (poles) in the domain.

    A pole is detected if:
    1. Q changes sign (crosses zero), OR
    2. min(|Q|) / max(|Q|) < relative_threshold (very small relative to maximum)
    """
    from rational_collocation import BernsteinBasis
    x_check = np.linspace(a, b, 200)
    Q_vals = np.array([BernsteinBasis.evaluate(result.Q_coeffs, xi, a, b) for xi in x_check])

    min_Q = np.min(Q_vals)
    max_Q = np.max(Q_vals)
    min_abs_Q = np.min(np.abs(Q_vals))
    max_abs_Q = np.max(np.abs(Q_vals))

    # Check for sign change
    if min_Q * max_Q < 0:
        return True, min_abs_Q

    # Check for near-zero relative to maximum
    if max_abs_Q > 0 and min_abs_Q / max_abs_Q < relative_threshold:
        return True, min_abs_Q

    return False, min_abs_Q


def main():
    """Run comprehensive benchmark suite"""

    # Test configurations
    # Use conservative degrees for initial demonstration
    # Note: Strong regularization keeps Q ≈ 1 (polynomial behavior)
    rational_degrees = [4, 6, 8]
    fd_grid_sizes = [10, 20, 40, 80, 160]

    # Get all test problems
    problems = [
        BVPBenchmarkSuite.smooth_poisson(),
        BVPBenchmarkSuite.smooth_polynomial_quartic(),
        BVPBenchmarkSuite.smooth_trig(),
        BVPBenchmarkSuite.smooth_exp(),
    ]

    all_results = []

    # Run benchmarks
    for problem in problems:
        results = benchmark_problem(problem, rational_degrees, fd_grid_sizes)
        all_results.extend(results)

    # Save results to JSON
    output_dir = Path(__file__).parent.parent / 'data'
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / 'bvp_rational_collocation_benchmark.json'

    # Convert to serializable format
    results_dict = {
        'metadata': {
            'description': 'BVP solver comparison: Rational collocation vs Polynomial FD',
            'rational_formulation': 'Quadratic form with explicit u and u\' unknowns',
            'solver': 'Levenberg-Marquardt for rational, direct solve for FD',
            'rational_degrees': rational_degrees,
            'fd_grid_sizes': fd_grid_sizes,
        },
        'results': [asdict(r) for r in all_results]
    }

    with open(output_file, 'w') as f:
        json.dump(results_dict, f, indent=2)

    print(f"\n{'='*70}")
    print(f"Benchmark complete!")
    print(f"Results saved to: {output_file}")
    print(f"Total tests: {len(all_results)}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
