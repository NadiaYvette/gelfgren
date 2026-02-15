#!/usr/bin/env python3
"""
Condensed Benchmark: Key problems and representative resolutions

Runs subset to demonstrate method comparison in reasonable time:
- 2 Helmholtz problems
- 2 Advection-Diffusion (including critical failure)
- 2 Reaction-Diffusion
- 2 Variable Coefficient (including critical failure)

Selected resolutions:
- FD: [20, 80]
- Spectral: [16, 32]
- Rational: [(6,3), (10,5)]
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / 'solvers'))
sys.path.insert(0, str(Path(__file__).parent / 'problems'))
sys.path.insert(0, str(Path(__file__).parent))

from run_extended_benchmark import (
    solve_with_polynomial_fd,
    solve_with_high_order_fd,
    solve_with_spectral,
    solve_with_rational,
    save_results,
    BenchmarkResult
)

from helmholtz_problems import helmholtz_sin_k1, helmholtz_exp_k4
from advection_diffusion_problems import advection_diffusion_mild, advection_diffusion_sharp
from reaction_diffusion_problems import reaction_diffusion_exponential_layers, reaction_diffusion_polynomial
from variable_coefficient_problems import variable_coefficient_polynomial, variable_coefficient_discontinuous


def run_condensed_benchmark():
    """Run condensed benchmark suite"""
    print("=" * 80)
    print("Condensed BVP Benchmark Suite")
    print("8 key problems × 6 methods × 2-3 resolutions")
    print("=" * 80)

    # Select representative problems
    problems = [
        helmholtz_sin_k1(),
        helmholtz_exp_k4(),
        advection_diffusion_mild(),
        advection_diffusion_sharp(),  # CRITICAL: expect rational failure
        reaction_diffusion_exponential_layers(),
        reaction_diffusion_polynomial(),
        variable_coefficient_polynomial(),
        variable_coefficient_discontinuous(),  # CRITICAL: expect rational failure
    ]

    # Condensed configuration
    fd_sizes = [20, 80]
    spectral_degrees = [16, 32]
    rational_degrees = [(6, 3), (10, 5)]

    all_results = []

    for i, problem in enumerate(problems, 1):
        print(f"\n{'=' * 80}")
        print(f"Problem {i}/{len(problems)}: {problem.name}")
        print(f"Type: {problem.operator_type}, Difficulty: {problem.expected_difficulty}")
        if "FAIL" in problem.description:
            print(f"*** {problem.description} ***")
        print(f"{'=' * 80}")
        print(f"{'Method':<40} {'DOF':<8} {'L2 Error':<15} {'Time (ms)':<12} {'Status':<10}")
        print("-" * 80)

        # Test each method
        for n_grid in fd_sizes:
            # 2nd-order FD
            r = solve_with_polynomial_fd(problem, n_grid)
            if r:
                all_results.append(r)
                print(f"{r.method_name + f' n={n_grid}':<40} {r.dof:<8} {r.l2_error:<15.2e} "
                      f"{r.solve_time_ms:<12.2f} {'OK' if r.success else 'FAIL'}")

            # 4th-order FD
            if n_grid >= 4:
                r = solve_with_high_order_fd(problem, n_grid, order=4)
                if r:
                    all_results.append(r)
                    print(f"{r.method_name + f' n={n_grid}':<40} {r.dof:<8} {r.l2_error:<15.2e} "
                          f"{r.solve_time_ms:<12.2f} {'OK' if r.success else 'FAIL'}")

        # 6th-order FD (only n=80)
        if 80 in fd_sizes:
            r = solve_with_high_order_fd(problem, 80, order=6)
            if r:
                all_results.append(r)
                print(f"{r.method_name + f' n=80':<40} {r.dof:<8} {r.l2_error:<15.2e} "
                      f"{r.solve_time_ms:<12.2f} {'OK' if r.success else 'FAIL'}")

        # Spectral methods
        for N in spectral_degrees:
            # Chebyshev
            r = solve_with_spectral(problem, N, 'chebyshev')
            if r:
                all_results.append(r)
                print(f"{'Chebyshev_Spectral' + f' N={N}':<40} {r.dof:<8} {r.l2_error:<15.2e} "
                      f"{r.solve_time_ms:<12.2f} {'OK' if r.success else 'FAIL'}")

            # Legendre
            r = solve_with_spectral(problem, N, 'legendre')
            if r:
                all_results.append(r)
                print(f"{'Legendre_Spectral' + f' N={N}':<40} {r.dof:<8} {r.l2_error:<15.2e} "
                      f"{r.solve_time_ms:<12.2f} {'OK' if r.success else 'FAIL'}")

        # Rational methods
        for n, m in rational_degrees:
            r = solve_with_rational(problem, n, m)
            if r:
                all_results.append(r)
                status = 'FAIL**' if not r.success or r.l2_error > 0.1 else 'OK'
                print(f"{f'Rational_Cleared[{n}/{m}]':<40} {r.dof:<8} {r.l2_error:<15.2e} "
                      f"{r.solve_time_ms:<12.2f} {status}")
                if r.notes:
                    print(f"  → {r.notes[:70]}")

    return all_results


def main():
    results = run_condensed_benchmark()

    # Save results
    output_dir = Path(__file__).parent.parent / 'data' / 'extended_benchmark'
    output_file = output_dir / 'condensed_benchmark_results.json'
    save_results(results, output_file)

    # Print summary
    print(f"\n{'=' * 80}")
    print("SUMMARY: Key Findings")
    print(f"{'=' * 80}\n")

    # Count failures
    rational_results = [r for r in results if 'Rational' in r.method_name]
    smooth_rationals = [r for r in rational_results if r.problem_name in [
        'Helmholtz_Sin_k1', 'Helmholtz_Exp_k4', 'ReacDiff_Exp_c10', 'ReacDiff_Polynomial_c4',
        'VarCoeff_Poly', 'AdvDiff_Mild_eps0.1'
    ]]
    nonsmooth_rationals = [r for r in rational_results if r.problem_name in [
        'AdvDiff_Sharp_eps0.001', 'VarCoeff_Discontinuous'
    ]]

    smooth_successes = sum(1 for r in smooth_rationals if r.success and r.l2_error < 1e-3)
    smooth_total = len(smooth_rationals)

    nonsmooth_failures = sum(1 for r in nonsmooth_rationals if not r.success or r.l2_error > 0.1)
    nonsmooth_total = len(nonsmooth_rationals)

    print(f"Rational Method Performance:")
    print(f"  Smooth problems: {smooth_successes}/{smooth_total} achieved < 1e-3 error")
    print(f"  Non-smooth problems: {nonsmooth_failures}/{nonsmooth_total} CATASTROPHIC FAILURES (>10% error)")
    print(f"\nConclusion: Rationals excel on smooth problems but catastrophically fail on")
    print(f"            boundary layers and discontinuous coefficients.")
    print(f"{'=' * 80}\n")


if __name__ == "__main__":
    main()
