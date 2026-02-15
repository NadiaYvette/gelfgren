#!/usr/bin/env python3
"""
Analysis of Extended Benchmark Results

Computes convergence rates, generates comparison plots, and identifies
best methods for each problem type.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict
from scipy import stats


def load_results(results_file):
    """Load benchmark results from JSON"""
    with open(results_file, 'r') as f:
        data = json.load(f)
    return data['results'], data.get('metadata', {})


def compute_convergence_rate(dofs, errors):
    """
    Compute convergence rate α where error ~ DOF^(-α).

    Uses log-log linear regression: log(error) = -α*log(DOF) + C
    """
    if len(dofs) < 2:
        return None

    # Filter out failures and zero errors
    valid = [(d, e) for d, e in zip(dofs, errors)
             if e > 0 and e < 1e10 and np.isfinite(e)]

    if len(valid) < 2:
        return None

    dofs_valid, errors_valid = zip(*valid)

    # Log-log linear regression
    log_dof = np.log(dofs_valid)
    log_err = np.log(errors_valid)

    slope, intercept, r_value, p_value, std_err = stats.linregress(log_dof, log_err)

    return -slope  # Convergence rate α


def analyze_method_on_problem(results, method_name, problem_name):
    """
    Analyze single method on single problem.

    Returns dict with convergence rate, errors, DOFs, etc.
    """
    # Filter results for this method and problem
    method_results = [r for r in results
                     if r['problem_name'] == problem_name
                     and method_name in r['method_name']]

    if not method_results:
        return None

    # Sort by DOF
    method_results.sort(key=lambda r: r['dof'])

    dofs = [r['dof'] for r in method_results]
    l2_errors = [r['l2_error'] for r in method_results]
    max_errors = [r['max_error'] for r in method_results]
    times = [r['solve_time_ms'] for r in method_results]
    successes = [r['success'] for r in method_results]

    # Compute convergence rate
    conv_rate = compute_convergence_rate(dofs, l2_errors)

    # Best result (minimum error among successful runs)
    successful_errors = [e for e, s in zip(l2_errors, successes) if s and e < 1e10]
    best_error = min(successful_errors) if successful_errors else float('inf')

    return {
        'method': method_name,
        'problem': problem_name,
        'dofs': dofs,
        'l2_errors': l2_errors,
        'max_errors': max_errors,
        'times': times,
        'successes': successes,
        'convergence_rate': conv_rate,
        'best_error': best_error,
        'num_tests': len(method_results)
    }


def generate_convergence_plot(results, problem_name, output_file):
    """
    Generate error vs DOF plot for all methods on one problem.
    """
    # Extract unique methods
    methods = list(set(r['method_name'].split('[')[0].split(' ')[0]
                      for r in results if r['problem_name'] == problem_name))

    fig, ax = plt.subplots(figsize=(10, 7))

    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    color_map = {
        'Polynomial_FD_2nd': colors[0],
        'Polynomial_FD_4th': colors[1],
        'Polynomial_FD_6th': colors[2],
        'Chebyshev_Spectral': colors[3],
        'Legendre_Spectral': colors[4],
        'Rational_Cleared': colors[5],
    }

    markers = {
        'Polynomial_FD_2nd': 'o',
        'Polynomial_FD_4th': 's',
        'Polynomial_FD_6th': '^',
        'Chebyshev_Spectral': 'D',
        'Legendre_Spectral': 'v',
        'Rational_Cleared': '*',
    }

    for method in sorted(set(r['method_name'] for r in results if r['problem_name'] == problem_name)):
        method_base = method.split('[')[0].split(' ')[0]

        analysis = analyze_method_on_problem(results, method_base, problem_name)
        if not analysis or not analysis['dofs']:
            continue

        dofs = analysis['dofs']
        errors = analysis['l2_errors']

        # Filter valid points
        valid = [(d, e) for d, e in zip(dofs, errors) if e > 0 and e < 1e10]
        if not valid:
            continue

        dofs_valid, errors_valid = zip(*valid)

        label = method_base.replace('_', ' ')
        if analysis['convergence_rate']:
            label += f" (α={analysis['convergence_rate']:.2f})"

        ax.loglog(dofs_valid, errors_valid,
                 marker=markers.get(method_base, 'o'),
                 color=color_map.get(method_base, 'gray'),
                 linewidth=2, markersize=8,
                 label=label, alpha=0.8)

    ax.set_xlabel('Degrees of Freedom (DOF)', fontsize=12)
    ax.set_ylabel('L² Error', fontsize=12)
    ax.set_title(f'{problem_name}: Error vs DOF', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=9, loc='best')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_file}")


def generate_convergence_rate_comparison(results, problem_name, output_file):
    """
    Generate bar chart comparing convergence rates across methods.
    """
    methods = ['Polynomial_FD_2nd', 'Polynomial_FD_4th', 'Polynomial_FD_6th',
               'Chebyshev_Spectral', 'Legendre_Spectral', 'Rational_Cleared']

    rates = []
    labels = []
    colors_list = []

    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    color_map = {
        'Polynomial_FD_2nd': colors[0],
        'Polynomial_FD_4th': colors[1],
        'Polynomial_FD_6th': colors[2],
        'Chebyshev_Spectral': colors[3],
        'Legendre_Spectral': colors[4],
        'Rational_Cleared': colors[5],
    }

    for method in methods:
        analysis = analyze_method_on_problem(results, method, problem_name)
        if analysis and analysis['convergence_rate']:
            rates.append(analysis['convergence_rate'])
            labels.append(method.replace('_', ' '))
            colors_list.append(color_map[method])

    if not rates:
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    bars = ax.bar(range(len(rates)), rates, color=colors_list, alpha=0.7, edgecolor='black')
    ax.set_xticks(range(len(rates)))
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylabel('Convergence Rate α', fontsize=12)
    ax.set_title(f'{problem_name}: Convergence Rates', fontsize=14, fontweight='bold')
    ax.grid(True, axis='y', alpha=0.3)

    # Add value labels on bars
    for bar, rate in zip(bars, rates):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
               f'{rate:.2f}',
               ha='center', va='bottom', fontsize=10)

    # Add reference lines
    ax.axhline(y=2, color='red', linestyle='--', alpha=0.5, label='2nd order')
    ax.axhline(y=4, color='orange', linestyle='--', alpha=0.5, label='4th order')
    ax.axhline(y=6, color='green', linestyle='--', alpha=0.5, label='6th order')
    ax.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_file}")


def generate_failure_comparison(results, output_file):
    """
    Generate visualization of method failures on critical problems.
    """
    critical_problems = [
        'AdvDiff_Sharp_eps0.001',
        'VarCoeff_Discontinuous'
    ]

    methods = ['Polynomial_FD_4th', 'Chebyshev_Spectral', 'Rational_Cleared']

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for idx, problem in enumerate(critical_problems):
        ax = axes[idx]

        method_errors = []
        method_labels = []
        colors_list = []

        color_map = {
            'Polynomial_FD_4th': 'steelblue',
            'Chebyshev_Spectral': 'orange',
            'Rational_Cleared': 'red'
        }

        for method in methods:
            analysis = analyze_method_on_problem(results, method, problem)
            if analysis and analysis['best_error'] < 1e10:
                method_errors.append(analysis['best_error'] * 100)  # Convert to percentage
                method_labels.append(method.replace('_', ' '))
                colors_list.append(color_map[method])

        if method_errors:
            bars = ax.bar(range(len(method_errors)), method_errors,
                         color=colors_list, alpha=0.7, edgecolor='black')
            ax.set_xticks(range(len(method_errors)))
            ax.set_xticklabels(method_labels, rotation=45, ha='right')
            ax.set_ylabel('Best L² Error (%)', fontsize=11)
            ax.set_title(problem.replace('_', ' '), fontsize=12, fontweight='bold')
            ax.set_yscale('log')
            ax.grid(True, axis='y', alpha=0.3)

            # Add 10% failure threshold line
            ax.axhline(y=10, color='darkred', linestyle='--', linewidth=2,
                      label='Failure Threshold (10%)')
            ax.legend(fontsize=9)

            # Annotate catastrophic failures
            for bar, err in zip(bars, method_errors):
                if err > 10:
                    ax.text(bar.get_x() + bar.get_width()/2., bar.get_height(),
                           'FAIL', ha='center', va='bottom',
                           fontsize=10, fontweight='bold', color='red')

    plt.suptitle('Critical Failure Cases: Rational vs Other Methods',
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_file}")


def generate_summary_table(results, output_file):
    """
    Generate summary table of best method for each problem.
    """
    problems = sorted(set(r['problem_name'] for r in results))
    methods = ['Polynomial_FD_2nd', 'Polynomial_FD_4th', 'Polynomial_FD_6th',
               'Chebyshev_Spectral', 'Legendre_Spectral', 'Rational_Cleared']

    with open(output_file, 'w') as f:
        f.write("# Extended Benchmark Analysis Summary\n\n")
        f.write("## Best Method for Each Problem\n\n")
        f.write("| Problem | Best Method | Error | Convergence Rate |\n")
        f.write("|---------|-------------|-------|------------------|\n")

        for problem in problems:
            best_method = None
            best_error = float('inf')
            best_rate = None

            for method in methods:
                analysis = analyze_method_on_problem(results, method, problem)
                if analysis and analysis['best_error'] < best_error:
                    best_error = analysis['best_error']
                    best_method = method
                    best_rate = analysis['convergence_rate']

            if best_method:
                rate_str = f"{best_rate:.2f}" if best_rate else "N/A"
                f.write(f"| {problem} | {best_method.replace('_', ' ')} | "
                       f"{best_error:.2e} | {rate_str} |\n")

        f.write("\n## Convergence Rate Summary\n\n")
        f.write("| Problem | 2nd FD | 4th FD | 6th FD | Cheby | Leg | Rational |\n")
        f.write("|---------|--------|--------|--------|-------|-----|----------|\n")

        for problem in problems:
            row = [problem]
            for method in methods:
                analysis = analyze_method_on_problem(results, method, problem)
                if analysis and analysis['convergence_rate']:
                    row.append(f"{analysis['convergence_rate']:.2f}")
                else:
                    row.append("N/A")
            f.write("| " + " | ".join(row) + " |\n")

    print(f"  Saved: {output_file}")


def main():
    """Run complete analysis"""
    print("=" * 80)
    print("Extended Benchmark Analysis")
    print("=" * 80)

    # Load results
    results_file = Path(__file__).parent.parent / 'data' / 'extended_benchmark' / 'condensed_benchmark_results.json'
    results, metadata = load_results(results_file)

    print(f"\nLoaded {len(results)} test results")

    # Create output directory
    output_dir = Path(__file__).parent.parent / 'reports' / 'figures' / 'extended_benchmark'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get unique problems
    problems = sorted(set(r['problem_name'] for r in results))
    print(f"Analyzing {len(problems)} problems")

    # Generate convergence plots for each problem
    print("\nGenerating convergence plots...")
    for problem in problems:
        output_file = output_dir / f'{problem}_convergence.pdf'
        generate_convergence_plot(results, problem, output_file)

    # Generate convergence rate comparisons
    print("\nGenerating convergence rate comparisons...")
    for problem in problems:
        output_file = output_dir / f'{problem}_rates.pdf'
        generate_convergence_rate_comparison(results, problem, output_file)

    # Generate failure comparison
    print("\nGenerating critical failure comparison...")
    failure_file = output_dir / 'critical_failures.pdf'
    generate_failure_comparison(results, failure_file)

    # Generate summary table
    print("\nGenerating summary table...")
    summary_file = output_dir / 'analysis_summary.md'
    generate_summary_table(results, summary_file)

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print(f"Figures saved to: {output_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
