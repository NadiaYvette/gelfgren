#!/usr/bin/env python3
"""Test critical failure cases: rationals on boundary layers and discontinuities"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / 'solvers'))
sys.path.insert(0, str(Path(__file__).parent / 'problems'))
sys.path.insert(0, str(Path(__file__).parent))

from advection_diffusion_problems import advection_diffusion_sharp
from variable_coefficient_problems import variable_coefficient_discontinuous
from run_extended_benchmark import solve_with_rational, solve_with_spectral

print("=" * 80)
print("CRITICAL FAILURE TESTS: Rationals on Non-Smooth Problems")
print("=" * 80)

# Test 1: Sharp boundary layer (ε = 0.001)
problem1 = advection_diffusion_sharp()
print(f"\nTest 1: {problem1.name}")
print(f"Description: {problem1.description}")
print(f"Expected: CATASTROPHIC FAILURE\n")

for n, m in [(6, 3), (10, 5), (16, 8)]:
    result = solve_with_rational(problem1, n, m)
    status = "**FAILED**" if not result.success or result.l2_error > 0.1 else "OK"
    print(f"Rational [{n}/{m}]: Error = {result.l2_error:.2e}, Success = {result.success} {status}")
    if result.notes:
        print(f"  Notes: {result.notes}")

# Compare with Chebyshev spectral (should also struggle but not as badly)
print("\nComparison with Chebyshev spectral:")
for N in [32, 64]:
    result = solve_with_spectral(problem1, N, 'chebyshev')
    if result:
        print(f"Chebyshev N={N}: Error = {result.l2_error:.2e}, Success = {result.success}")

# Test 2: Discontinuous coefficient
print(f"\n{'=' * 80}")
problem2 = variable_coefficient_discontinuous()
print(f"\nTest 2: {problem2.name}")
print(f"Description: {problem2.description}")
print(f"Expected: CATASTROPHIC FAILURE\n")

for n, m in [(6, 3), (10, 5), (16, 8)]:
    result = solve_with_rational(problem2, n, m)
    status = "**FAILED**" if not result.success or result.l2_error > 0.1 else "OK"
    print(f"Rational [{n}/{m}]: Error = {result.l2_error:.2e}, Success = {result.success} {status}")
    if result.notes:
        print(f"  Notes: {result.notes}")

print(f"\n{'=' * 80}")
print("Summary: These critical tests demonstrate rational method limitations")
print("=" * 80)
