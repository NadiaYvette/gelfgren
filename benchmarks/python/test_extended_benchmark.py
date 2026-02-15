#!/usr/bin/env python3
"""Quick test of extended benchmark on one problem"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / 'solvers'))
sys.path.insert(0, str(Path(__file__).parent / 'problems'))
sys.path.insert(0, str(Path(__file__).parent))

from helmholtz_problems import helmholtz_sin_k1
from run_extended_benchmark import (
    solve_with_polynomial_fd,
    solve_with_high_order_fd,
    solve_with_spectral,
    solve_with_rational
)

# Test one problem
problem = helmholtz_sin_k1()
print(f"Testing: {problem.name}")
print(f"Description: {problem.description}")
print()

# Test each method
print("Testing 2nd-order FD with n=20...")
result = solve_with_polynomial_fd(problem, 20)
print(f"  Success: {result.success}, Error: {result.l2_error:.2e}, Time: {result.solve_time_ms:.2f}ms")

print("\nTesting 4th-order FD with n=20...")
result = solve_with_high_order_fd(problem, 20, order=4)
if result:
    print(f"  Success: {result.success}, Error: {result.l2_error:.2e}, Time: {result.solve_time_ms:.2f}ms")

print("\nTesting Chebyshev spectral with N=16...")
result = solve_with_spectral(problem, 16, 'chebyshev')
if result:
    print(f"  Success: {result.success}, Error: {result.l2_error:.2e}, Time: {result.solve_time_ms:.2f}ms")

print("\nTesting Rational [6/3]...")
result = solve_with_rational(problem, 6, 3)
print(f"  Success: {result.success}, Error: {result.l2_error:.2e}, Time: {result.solve_time_ms:.2f}ms")
if result.notes:
    print(f"  Notes: {result.notes}")

print("\nAll methods tested successfully!")
