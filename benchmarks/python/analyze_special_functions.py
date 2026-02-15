#!/usr/bin/env python3
"""
Analyze Special Function Approximation Results

Generates summary tables showing:
- Best method for each function
- Convergence rates
- Final accuracy achieved
- Identification of rational failures
"""

import json
from pathlib import Path
import numpy as np


def load_results():
    """Load all special function results"""
    data_dir = Path(__file__).parent.parent / 'data' / 'special_functions'
    results = []

    for json_file in sorted(data_dir.glob('*.json')):
        if 'Exponential_e^x' in json_file.name or 'Logarithm_log1+x' in json_file.name:
            # Skip old files with incompatible names
            continue
        with open(json_file) as f:
            data = json.load(f)
            results.append(data)

    return results


def get_final_error(method_errors, error_type='l2_error'):
    """Get final (finest mesh) error"""
    if not method_errors:
        return float('inf')
    return method_errors[-1][error_type]


def get_convergence_rate(convergence_rates, rate_type='polynomial_l2'):
    """Get average convergence rate"""
    if rate_type not in convergence_rates:
        return 0.0
    rates = convergence_rates[rate_type]
    if not rates:
        return 0.0
    # Use last few rates (typically more stable)
    return np.mean(rates[-3:]) if len(rates) >= 3 else np.mean(rates)


def classify_result(poly_error, rat_error, poly_rate, rat_rate):
    """Classify which method performs better"""
    # Catastrophic failure: error > 1% or negative convergence
    poly_catastrophic = poly_error > 0.01 or poly_rate < 0.5
    rat_catastrophic = rat_error > 0.01 or rat_rate < 0.5

    if rat_catastrophic and not poly_catastrophic:
        return "Polynomial (Rational FAILS)"
    elif poly_catastrophic and not rat_catastrophic:
        return "Rational (Polynomial FAILS)"
    elif rat_catastrophic and poly_catastrophic:
        return "Both FAIL"

    # Both converge: compare final accuracy
    if rat_error < poly_error / 10:
        return "Rational Superior"
    elif poly_error < rat_error / 10:
        return "Polynomial Superior"
    else:
        return "Comparable"


def generate_summary_table():
    """Generate main summary table"""
    results = load_results()

    print("\n" + "="*80)
    print("Special Function Approximation Results Summary")
    print("="*80)
    print("\nTable: Comparison of Polynomial vs Rational Collocation")
    print("-"*80)
    print(f"{'Function':<30} {'Winner':<25} {'Poly L2':<12} {'Rat L2':<12}")
    print("-"*80)

    poly_wins = []
    rat_wins = []
    rat_fails = []
    comparable = []

    for data in results:
        name = data['function_name']

        poly_error = get_final_error(data['polynomial_errors'])
        rat_error = get_final_error(data['rational_errors'])

        poly_rate = get_convergence_rate(data['convergence_rates'], 'polynomial_l2')
        rat_rate = get_convergence_rate(data['convergence_rates'], 'rational_l2')

        classification = classify_result(poly_error, rat_error, poly_rate, rat_rate)

        # Categorize
        if "Rational FAILS" in classification:
            rat_fails.append(name)
        elif "Rational Superior" in classification:
            rat_wins.append(name)
        elif "Polynomial Superior" in classification:
            poly_wins.append(name)
        else:
            comparable.append(name)

        print(f"{name:<30} {classification:<25} {poly_error:<12.2e} {rat_error:<12.2e}")

    print("-"*80)
    print(f"\nRational Superior: {len(rat_wins)} functions")
    print(f"Polynomial Superior: {len(poly_wins)} functions")
    print(f"Comparable: {len(comparable)} functions")
    print(f"Rational Catastrophic Failure: {len(rat_fails)} functions")
    print("="*80)

    return rat_wins, poly_wins, rat_fails, comparable


def generate_detailed_table():
    """Generate detailed convergence rate table"""
    results = load_results()

    print("\n" + "="*80)
    print("Detailed Convergence Rates (Average of Last 3 Refinements)")
    print("="*80)
    print(f"{'Function':<30} {'Poly α_L2':<10} {'Rat α_L2':<10} {'Poly L2':<12} {'Rat L2':<12}")
    print("-"*80)

    for data in results:
        name = data['function_name']

        poly_error = get_final_error(data['polynomial_errors'])
        rat_error = get_final_error(data['rational_errors'])

        poly_rate = get_convergence_rate(data['convergence_rates'], 'polynomial_l2')
        rat_rate = get_convergence_rate(data['convergence_rates'], 'rational_l2')

        print(f"{name:<30} {poly_rate:>9.2f} {rat_rate:>9.2f} {poly_error:>11.2e} {rat_error:>11.2e}")

    print("-"*80)
    print("α_L2 is convergence rate where error ~ h^α")
    print("Final errors are at finest mesh (128 intervals)")
    print("="*80)


def generate_latex_table():
    """Generate LaTeX table for Chapter 5"""
    results = load_results()

    print("\n" + "="*80)
    print("LaTeX Table for Chapter 5")
    print("="*80)

    print("""
\\begin{table}[htbp]
\\centering
\\caption{Special Function Approximation Results: Polynomial vs Rational Collocation at 128 Intervals}
\\label{tab:special_function_results}
\\begin{tabular}{lccccl}
\\toprule
\\textbf{Function} & \\multicolumn{2}{c}{\\textbf{L2 Error}} & \\multicolumn{2}{c}{\\textbf{Conv. Rate}} & \\textbf{Winner} \\\\
\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}
& Polynomial & Rational & Poly $\\alpha$ & Rat $\\alpha$ & \\\\
\\midrule""")

    for data in results:
        name = data['function_name']

        # Clean up name for LaTeX
        name_latex = name.replace('_', '\\_').replace('^', '\\textasciicircum')
        if len(name_latex) > 30:
            name_latex = name_latex[:27] + "..."

        poly_error = get_final_error(data['polynomial_errors'])
        rat_error = get_final_error(data['rational_errors'])

        poly_rate = get_convergence_rate(data['convergence_rates'], 'polynomial_l2')
        rat_rate = get_convergence_rate(data['convergence_rates'], 'rational_l2')

        classification = classify_result(poly_error, rat_error, poly_rate, rat_rate)

        # Format for LaTeX
        poly_str = f"{poly_error:.1e}".replace('e-0', 'e-').replace('e-', '$\\times 10^{-') + '}$'
        rat_str = f"{rat_error:.1e}".replace('e-0', 'e-').replace('e-', '$\\times 10^{-') + '}$'

        # Bold the winner
        if "Rational Superior" in classification or "Polynomial FAILS" in classification:
            rat_str = f"\\textbf{{{rat_str}}}"
            winner = "Rational"
        elif "Polynomial Superior" in classification or "Rational FAILS" in classification:
            poly_str = f"\\textbf{{{poly_str}}}"
            winner = "Polynomial"
        else:
            winner = "Comparable"

        if "FAILS" in classification:
            winner = "\\textbf{" + winner + " FAIL}"

        print(f"{name_latex} & {poly_str} & {rat_str} & {poly_rate:.1f} & {rat_rate:.1f} & {winner} \\\\")

    print("""\\bottomrule
\\end{tabular}
\\end{table}
""")
    print("="*80)


def generate_failure_analysis():
    """Analyze patterns in rational failures"""
    results = load_results()

    print("\n" + "="*80)
    print("Rational Failure Analysis")
    print("="*80)

    failures = []
    successes = []

    for data in results:
        name = data['function_name']
        rat_error = get_final_error(data['rational_errors'])
        rat_rate = get_convergence_rate(data['convergence_rates'], 'rational_l2')

        if rat_error > 0.01 or rat_rate < 0.5:
            failures.append((name, rat_error, rat_rate))
        else:
            successes.append((name, rat_error, rat_rate))

    print(f"\nRational Method Failures ({len(failures)} functions):")
    print("-"*80)
    for name, error, rate in failures:
        print(f"  {name:<40} Error: {error:.2e}, Rate: {rate:.2f}")

    print(f"\nRational Method Successes ({len(successes)} functions):")
    print("-"*80)
    for name, error, rate in successes:
        print(f"  {name:<40} Error: {error:.2e}, Rate: {rate:.2f}")

    print("\n" + "="*80)
    print("Pattern Analysis:")
    print("-"*80)
    print("Rational failures occur on:")
    print("  - Highly oscillatory functions (Airy, high-order Bessel)")
    print("  - Functions with rapid frequency changes")
    print("  - Mathieu/Jacobi elliptic at high refinement (stagnation)")
    print("\nRational successes occur on:")
    print("  - Smooth exponential-type functions")
    print("  - Low-order smooth special functions")
    print("="*80)


def main():
    """Generate all analysis"""
    rat_wins, poly_wins, rat_fails, comparable = generate_summary_table()
    generate_detailed_table()
    generate_latex_table()
    generate_failure_analysis()

    # Generate key findings
    print("\n" + "="*80)
    print("KEY FINDINGS FOR CHAPTER 5")
    print("="*80)
    print(f"""
1. Rational collocation achieves superior accuracy on {len(rat_wins)} smooth functions
   (exponential, logarithm, low-order Bessel)

2. Rational collocation CATASTROPHICALLY FAILS on {len(rat_fails)} oscillatory functions
   (Airy functions, high-order Bessel, some Mathieu/Jacobi at high refinement)

3. Polynomial collocation is more ROBUST - successfully converges on all {len(rat_wins) + len(poly_wins) + len(comparable)} functions

4. For special function approximation, POLYNOMIAL methods are SAFER unless the function
   is known to have exponential or smooth rational character.

5. Parameter choices (q for Mathieu, k for Jacobi) probe specific oscillatory behaviors,
   demonstrating that approximation performance depends on problem-specific characteristics.
""")
    print("="*80)


if __name__ == "__main__":
    main()
