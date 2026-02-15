#!/usr/bin/env python3
"""
Extended Benchmark Suite: Comprehensive BVP Method Comparison

Compares 6 methods across 20+ test problems:
- Polynomial FD: 2nd, 4th, 6th order
- Spectral: Chebyshev, Legendre
- Rational: Cleared form collocation

Problem types:
- Helmholtz (5 tests)
- Advection-Diffusion (5 tests) - includes CRITICAL failure tests
- Reaction-Diffusion (5 tests)
- Variable Coefficient (5 tests) - includes CRITICAL failure tests
"""

import sys
from pathlib import Path

# Add paths for imports
sys.path.insert(0, str(Path(__file__).parent / 'solvers'))
sys.path.insert(0, str(Path(__file__).parent / 'problems'))
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import json
import time
from dataclasses import dataclass, asdict
from typing import List, Dict

# Import solvers
from spectral_collocation import ChebyshevSpectralSolver, LegendreSpectralSolver
from finite_differences import FourthOrderFD, SixthOrderFD
from rational_collocation_cleared import RationalCollocationCleared, QConstraintType
from compare_bvp_methods import PolynomialFiniteDifference

# Import problems
from helmholtz_problems import get_all_helmholtz_problems
from advection_diffusion_problems import get_all_advection_diffusion_problems
from reaction_diffusion_problems import get_all_reaction_diffusion_problems
from variable_coefficient_problems import get_all_variable_coefficient_problems


@dataclass
class BenchmarkResult:
    """Results from a single method on a single problem"""
    problem_name: str
    problem_type: str
    method_name: str
    degree_or_points: int  # N for spectral, n for FD, n for rational
    dof: int  # Total degrees of freedom
    max_error: float
    l2_error: float
    solve_time_ms: float
    success: bool
    residual_norm: float = 0.0
    notes: str = ""


def solve_with_polynomial_fd(problem, n_grid):
    """Solve using 2nd-order polynomial finite difference"""
    solver = PolynomialFiniteDifference(
        f=problem.forcing,
        a=problem.domain[0],
        b=problem.domain[1],
        alpha=problem.boundary_conditions['alpha'],
        beta=problem.boundary_conditions['beta']
    )

    start_time = time.time()
    try:
        x, u, success = solver.solve(n_grid)
        solve_time = time.time() - start_time

        # Evaluate error
        x_eval = np.linspace(problem.domain[0], problem.domain[1], 1000)
        u_computed = np.interp(x_eval, x, u)
        u_exact = problem.exact(x_eval)
        error = np.abs(u_computed - u_exact)

        return BenchmarkResult(
            problem_name=problem.name,
            problem_type=problem.operator_type,
            method_name='Polynomial_FD_2nd',
            degree_or_points=n_grid,
            dof=n_grid + 2,
            max_error=np.max(error),
            l2_error=np.sqrt(np.mean(error**2)),
            solve_time_ms=solve_time * 1000,
            success=success
        )
    except Exception as e:
        return BenchmarkResult(
            problem_name=problem.name,
            problem_type=problem.operator_type,
            method_name='Polynomial_FD_2nd',
            degree_or_points=n_grid,
            dof=n_grid + 2,
            max_error=float('inf'),
            l2_error=float('inf'),
            solve_time_ms=0.0,
            success=False,
            notes=f"Error: {str(e)[:50]}"
        )


def solve_with_high_order_fd(problem, n_grid, order=4):
    """Solve using 4th or 6th order finite difference"""
    # Get coefficients from problem
    coeffs = problem.coeffs
    epsilon = coeffs.get('epsilon', 1.0)
    b_coeff = coeffs.get('b_coeff', 0.0)
    c_coeff = coeffs.get('c_coeff', 0.0)
    k_squared = coeffs.get('k_squared', 0.0)
    a_func = coeffs.get('a_func', None)

    # Skip variable coefficient for now (not implemented in FD solvers)
    if a_func is not None:
        return None

    SolverClass = FourthOrderFD if order == 4 else SixthOrderFD
    method_name = f'Polynomial_FD_{order}th'

    solver = SolverClass(
        f=problem.forcing,
        a=problem.domain[0],
        b=problem.domain[1],
        alpha=problem.boundary_conditions['alpha'],
        beta=problem.boundary_conditions['beta'],
        epsilon=epsilon,
        b_coeff=b_coeff,
        c_coeff=c_coeff,
        k_squared=k_squared
    )

    start_time = time.time()
    try:
        result = solver.solve(n_grid)
        solve_time = time.time() - start_time

        # Evaluate error
        x_eval = np.linspace(problem.domain[0], problem.domain[1], 1000)
        u_computed = np.interp(x_eval, result.x, result.u)
        u_exact = problem.exact(x_eval)
        error = np.abs(u_computed - u_exact)

        return BenchmarkResult(
            problem_name=problem.name,
            problem_type=problem.operator_type,
            method_name=method_name,
            degree_or_points=n_grid,
            dof=n_grid + 2,
            max_error=np.max(error),
            l2_error=np.sqrt(np.mean(error**2)),
            solve_time_ms=solve_time * 1000,
            success=result.success,
            residual_norm=result.residual_norm
        )
    except Exception as e:
        return BenchmarkResult(
            problem_name=problem.name,
            problem_type=problem.operator_type,
            method_name=method_name,
            degree_or_points=n_grid,
            dof=n_grid + 2,
            max_error=float('inf'),
            l2_error=float('inf'),
            solve_time_ms=0.0,
            success=False,
            notes=f"Error: {str(e)[:50]}"
        )


def solve_with_spectral(problem, N, spectral_type='chebyshev'):
    """Solve using Chebyshev or Legendre spectral collocation"""
    # Get coefficients from problem
    coeffs = problem.coeffs
    epsilon = coeffs.get('epsilon', 1.0)
    b_coeff = coeffs.get('b_coeff', 0.0)
    c_coeff = coeffs.get('c_coeff', 0.0)
    k_squared = coeffs.get('k_squared', 0.0)
    a_func = coeffs.get('a_func', None)

    # Skip variable coefficient for now
    if a_func is not None:
        return None

    SolverClass = ChebyshevSpectralSolver if spectral_type == 'chebyshev' else LegendreSpectralSolver
    method_name = f'{spectral_type.capitalize()}_Spectral'

    solver = SolverClass(
        f=problem.forcing,
        a=problem.domain[0],
        b=problem.domain[1],
        alpha=problem.boundary_conditions['alpha'],
        beta=problem.boundary_conditions['beta'],
        epsilon=epsilon,
        b_coeff=b_coeff,
        c_coeff=c_coeff,
        k_squared=k_squared
    )

    start_time = time.time()
    try:
        result = solver.solve(N)
        solve_time = time.time() - start_time

        # Evaluate error
        x_eval = np.linspace(problem.domain[0], problem.domain[1], 1000)
        u_computed = np.interp(x_eval, result.x, result.u)
        u_exact = problem.exact(x_eval)
        error = np.abs(u_computed - u_exact)

        return BenchmarkResult(
            problem_name=problem.name,
            problem_type=problem.operator_type,
            method_name=method_name,
            degree_or_points=N,
            dof=N + 1,
            max_error=np.max(error),
            l2_error=np.sqrt(np.mean(error**2)),
            solve_time_ms=solve_time * 1000,
            success=result.success,
            residual_norm=result.residual_norm
        )
    except Exception as e:
        return BenchmarkResult(
            problem_name=problem.name,
            problem_type=problem.operator_type,
            method_name=method_name,
            degree_or_points=N,
            dof=N + 1,
            max_error=float('inf'),
            l2_error=float('inf'),
            solve_time_ms=0.0,
            success=False,
            notes=f"Error: {str(e)[:50]}"
        )


def solve_with_rational(problem, n, m):
    """Solve using rational collocation cleared form"""
    # Get coefficients from problem
    coeffs = problem.coeffs
    epsilon = coeffs.get('epsilon', 1.0)
    b_coeff = coeffs.get('b_coeff', 0.0)
    c_coeff = coeffs.get('c_coeff', 0.0)
    k_squared = coeffs.get('k_squared', 0.0)
    a_func = coeffs.get('a_func', None)

    solver = RationalCollocationCleared(
        f=problem.forcing,
        a=problem.domain[0],
        b=problem.domain[1],
        alpha=problem.boundary_conditions['alpha'],
        beta=problem.boundary_conditions['beta'],
        n=n,
        m=m,
        q_constraint=QConstraintType.ENDPOINT,
        epsilon=epsilon,
        b_coeff=b_coeff,
        c_coeff=c_coeff,
        k_squared=k_squared,
        a_func=a_func
    )

    start_time = time.time()
    try:
        result = solver.solve(method='trf', verbose=False)
        solve_time = time.time() - start_time

        # Evaluate error
        x_eval = np.linspace(problem.domain[0], problem.domain[1], 1000)
        u_computed = solver.evaluate_solution(result, x_eval)
        u_exact = problem.exact(x_eval)
        error = np.abs(u_computed - u_exact)

        return BenchmarkResult(
            problem_name=problem.name,
            problem_type=problem.operator_type,
            method_name=f'Rational_Cleared[{n}/{m}]',
            degree_or_points=n,
            dof=n + m,
            max_error=np.max(error),
            l2_error=np.sqrt(np.mean(error**2)),
            solve_time_ms=solve_time * 1000,
            success=result.success,
            residual_norm=result.residual_norm
        )
    except Exception as e:
        return BenchmarkResult(
            problem_name=problem.name,
            problem_type=problem.operator_type,
            method_name=f'Rational_Cleared[{n}/{m}]',
            degree_or_points=n,
            dof=n + m,
            max_error=float('inf'),
            l2_error=float('inf'),
            solve_time_ms=0.0,
            success=False,
            notes=f"Error: {str(e)[:100]}"
        )


def run_benchmark_suite():
    """Run comprehensive benchmark suite"""
    print("=" * 80)
    print("Extended BVP Benchmark Suite")
    print("=" * 80)

    # Collect all problems
    all_problems = []
    all_problems.extend(get_all_helmholtz_problems())
    all_problems.extend(get_all_advection_diffusion_problems())
    all_problems.extend(get_all_reaction_diffusion_problems())
    all_problems.extend(get_all_variable_coefficient_problems())

    print(f"\nTotal problems: {len(all_problems)}")
    print(f"  Helmholtz: 5")
    print(f"  Advection-Diffusion: 5")
    print(f"  Reaction-Diffusion: 5")
    print(f"  Variable Coefficient: 5")

    # Configuration
    fd_grid_sizes = [10, 20, 40, 80, 160]
    spectral_degrees = [8, 16, 32, 64]
    rational_degrees = [(4, 2), (6, 3), (8, 4), (10, 5)]

    all_results = []

    for i, problem in enumerate(all_problems, 1):
        print(f"\n{'=' * 80}")
        print(f"Problem {i}/{len(all_problems)}: {problem.name}")
        print(f"Type: {problem.operator_type}, Difficulty: {problem.expected_difficulty}")
        print(f"{'=' * 80}")
        print(f"{'Method':<35} {'DOF':<8} {'Max Error':<15} {'L2 Error':<15} {'Time (ms)':<12}")
        print("-" * 80)

        # Test 2nd-order FD
        for n_grid in fd_grid_sizes:
            result = solve_with_polynomial_fd(problem, n_grid)
            if result:
                all_results.append(result)
                print(f"{result.method_name:<35} {result.dof:<8} {result.max_error:<15.2e} "
                      f"{result.l2_error:<15.2e} {result.solve_time_ms:<12.2f}")

        # Test 4th-order FD
        for n_grid in fd_grid_sizes:
            if n_grid >= 4:  # Need at least 4 interior points
                result = solve_with_high_order_fd(problem, n_grid, order=4)
                if result:
                    all_results.append(result)
                    print(f"{result.method_name:<35} {result.dof:<8} {result.max_error:<15.2e} "
                          f"{result.l2_error:<15.2e} {result.solve_time_ms:<12.2f}")

        # Test 6th-order FD
        for n_grid in [40, 80, 160]:  # Need at least 6 interior points
            if n_grid >= 6:
                result = solve_with_high_order_fd(problem, n_grid, order=6)
                if result:
                    all_results.append(result)
                    print(f"{result.method_name:<35} {result.dof:<8} {result.max_error:<15.2e} "
                          f"{result.l2_error:<15.2e} {result.solve_time_ms:<12.2f}")

        # Test Chebyshev spectral
        for N in spectral_degrees:
            result = solve_with_spectral(problem, N, 'chebyshev')
            if result:
                all_results.append(result)
                print(f"{result.method_name:<35} {result.dof:<8} {result.max_error:<15.2e} "
                      f"{result.l2_error:<15.2e} {result.solve_time_ms:<12.2f}")

        # Test Legendre spectral
        for N in spectral_degrees:
            result = solve_with_spectral(problem, N, 'legendre')
            if result:
                all_results.append(result)
                print(f"{result.method_name:<35} {result.dof:<8} {result.max_error:<15.2e} "
                      f"{result.l2_error:<15.2e} {result.solve_time_ms:<12.2f}")

        # Test rational collocation
        for n, m in rational_degrees:
            result = solve_with_rational(problem, n, m)
            if result:
                all_results.append(result)
                status = "FAILED" if not result.success or result.max_error > 1e-2 else ""
                print(f"{result.method_name:<35} {result.dof:<8} {result.max_error:<15.2e} "
                      f"{result.l2_error:<15.2e} {result.solve_time_ms:<12.2f} {status}")
                if result.notes:
                    print(f"  Note: {result.notes}")

    return all_results


def save_results(results: List[BenchmarkResult], output_file: Path):
    """Save results to JSON"""
    results_dict = {
        'metadata': {
            'description': 'Extended benchmark: 6 methods × 20 problems',
            'methods': [
                'Polynomial_FD_2nd',
                'Polynomial_FD_4th',
                'Polynomial_FD_6th',
                'Chebyshev_Spectral',
                'Legendre_Spectral',
                'Rational_Cleared'
            ],
            'problem_types': ['helmholtz', 'advection_diffusion', 'reaction_diffusion', 'variable_coefficient']
        },
        'results': [asdict(r) for r in results]
    }

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(results_dict, f, indent=2)

    print(f"\n{'=' * 80}")
    print(f"Results saved to: {output_file}")
    print(f"Total tests: {len(results)}")
    print(f"{'=' * 80}")


def main():
    results = run_benchmark_suite()

    # Save results
    output_dir = Path(__file__).parent.parent / 'data' / 'extended_benchmark'
    output_file = output_dir / 'extended_benchmark_results.json'
    save_results(results, output_file)


if __name__ == "__main__":
    main()
