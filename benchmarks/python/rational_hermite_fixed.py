#!/usr/bin/env python3
"""
Robust Rational Hermite Interpolation

Constructs [n/m] rational P(x)/Q(x) by matching function and derivatives
at interval endpoints using least-squares optimization.

For [n/m] rational matching k derivatives at each of 2 points:
- Conditions: 2(k+1)
- Unknowns: (n+1) + (m+1) - 1 = n+m+1 (after normalization)
- Require: 2(k+1) = n+m+1, giving k = (n+m-1)/2

Examples:
- [3/2]: k=2 (match u, u', u'' at both ends)
- [5/4]: k=4 (match u→u'''' at both ends)
- [7/6]: k=6 (match u→u⁶ at both ends)
"""

import numpy as np
from scipy.optimize import minimize, least_squares
from typing import Callable, List, Tuple


def compute_derivatives(func: Callable, x: float, max_order: int) -> np.ndarray:
    """Compute function and derivatives using finite differences"""
    h = 1e-6
    derivs = [func(x)]

    if max_order >= 1:
        derivs.append((func(x + h) - func(x - h)) / (2 * h))
    if max_order >= 2:
        derivs.append((func(x + h) - 2*func(x) + func(x - h)) / h**2)
    if max_order >= 3:
        d3 = (-func(x + 2*h) + 2*func(x + h) - 2*func(x - h) + func(x - 2*h)) / (2 * h**3)
        derivs.append(d3)
    if max_order >= 4:
        d4 = (func(x + 2*h) - 4*func(x + h) + 6*func(x) - 4*func(x - h) + func(x - 2*h)) / h**4
        derivs.append(d4)
    if max_order >= 5:
        d5 = (func(x + 3*h) - 6*func(x + 2*h) + 15*func(x + h) - 20*func(x) +
              15*func(x - h) - 6*func(x - 2*h) + func(x - 3*h)) / (2 * h**5)
        derivs.append(d5)
    if max_order >= 6:
        d6 = (func(x + 3*h) - 9*func(x + 2*h) + 45*func(x + h) - 45*func(x - h) +
              9*func(x - 2*h) - func(x - 3*h)) / (60 * h**6)
        derivs.append(d6)

    return np.array(derivs)


def eval_rational(p_coeffs: np.ndarray, q_coeffs: np.ndarray, x: float) -> float:
    """Evaluate rational P(x)/Q(x)"""
    P = sum(p_coeffs[i] * x**i for i in range(len(p_coeffs)))
    Q = sum(q_coeffs[i] * x**i for i in range(len(q_coeffs)))
    return P / Q if abs(Q) > 1e-12 else np.nan


def eval_rational_deriv(p_coeffs: np.ndarray, q_coeffs: np.ndarray, x: float, order: int) -> float:
    """
    Evaluate k-th derivative of rational P(x)/Q(x) using Faà di Bruno formula.

    For order 1: (P/Q)' = (P'Q - PQ')/Q²
    For higher orders, use numerical differentiation of lower order.
    """
    if order == 0:
        return eval_rational(p_coeffs, q_coeffs, x)

    # Use numerical differentiation for derivatives
    h = 1e-7
    if order == 1:
        return (eval_rational(p_coeffs, q_coeffs, x + h) -
                eval_rational(p_coeffs, q_coeffs, x - h)) / (2 * h)
    else:
        return (eval_rational_deriv(p_coeffs, q_coeffs, x + h, order - 1) -
                eval_rational_deriv(p_coeffs, q_coeffs, x - h, order - 1)) / (2 * h)


def construct_rational_hermite_robust(a: float, b: float,
                                       values_a: np.ndarray,
                                       values_b: np.ndarray,
                                       deg_num: int, deg_den: int,
                                       use_interior_pole_penalty: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """
    Construct [deg_num/deg_den] rational matching Hermite data.

    Uses least-squares with pole penalties and multiple initializations for robustness.
    """
    K = len(values_a) - 1  # Number of derivatives matched
    n_unknowns = deg_num + deg_den + 1  # After Q normalization

    # Strategy: Represent as P(x) = sum p_i x^i, Q(x) = 1 + sum q_i x^i
    # This automatically enforces Q(a) ≠ 0 (approximately)

    def unpack_coeffs(x):
        """Unpack optimization variables into P and Q coefficients"""
        p_coeffs = x[:deg_num + 1]
        q_coeffs = np.ones(deg_den + 1)
        if deg_den > 0:
            q_coeffs[1:] = x[deg_num + 1:]
        return p_coeffs, q_coeffs

    def residuals(x):
        """Compute residuals for Hermite interpolation conditions"""
        p_coeffs, q_coeffs = unpack_coeffs(x)
        res = []

        # Conditions at a
        for k in range(len(values_a)):
            try:
                r_val = eval_rational_deriv(p_coeffs, q_coeffs, a, k)
                if np.isnan(r_val):
                    res.append(1e10)
                else:
                    res.append(r_val - values_a[k])
            except:
                res.append(1e10)

        # Conditions at b
        for k in range(len(values_b)):
            try:
                r_val = eval_rational_deriv(p_coeffs, q_coeffs, b, k)
                if np.isnan(r_val):
                    res.append(1e10)
                else:
                    res.append(r_val - values_b[k])
            except:
                res.append(1e10)

        # Pole penalty: penalize Q ≈ 0 in interior
        if use_interior_pole_penalty:
            x_interior = np.linspace(a, b, 7)[1:-1]  # 5 interior points
            for xi in x_interior:
                Q_val = sum(q_coeffs[i] * xi**i for i in range(len(q_coeffs)))
                if abs(Q_val) < 0.1:
                    res.append(100 * (0.1 - abs(Q_val)))

        return np.array(res)

    # Initial guess 1: Polynomial interpolant with Q = 1
    x0_poly = np.zeros(n_unknowns)
    try:
        # Simple polynomial fit
        p_fit = np.polyfit([a, b], [values_a[0], values_b[0]], min(deg_num, 1))
        for i, coeff in enumerate(reversed(p_fit)):
            if i < deg_num + 1:
                x0_poly[i] = coeff
    except:
        x0_poly[:deg_num + 1] = [values_a[0]] + [0] * deg_num

    # Try least squares
    try:
        result = least_squares(residuals, x0_poly, method='lm', ftol=1e-12, xtol=1e-12, max_nfev=20000)
        if result.cost < 1e-6:  # Good solution
            p_coeffs, q_coeffs = unpack_coeffs(result.x)
            return p_coeffs, q_coeffs
    except:
        pass

    # Initial guess 2: Small random perturbation
    x0_rand = x0_poly + np.random.randn(n_unknowns) * 0.01
    try:
        result = least_squares(residuals, x0_rand, method='lm', ftol=1e-12, xtol=1e-12, max_nfev=20000)
        if result.cost < 1e-6:
            p_coeffs, q_coeffs = unpack_coeffs(result.x)
            return p_coeffs, q_coeffs
    except:
        pass

    # Initial guess 3: Different scale
    x0_scaled = x0_poly.copy()
    x0_scaled[deg_num + 1:] = 0.1  # Small Q perturbations
    try:
        result = least_squares(residuals, x0_scaled, method='lm', ftol=1e-10, xtol=1e-10, max_nfev=20000)
        if result.cost < 1e-4:
            p_coeffs, q_coeffs = unpack_coeffs(result.x)
            return p_coeffs, q_coeffs
    except:
        pass

    # Fallback: Return polynomial approximation
    print(f"Warning: Rational fit failed on [{a:.3f}, {b:.3f}], using polynomial fallback")
    p_coeffs = x0_poly[:deg_num + 1]
    q_coeffs = np.ones(deg_den + 1)
    return p_coeffs, q_coeffs


class PiecewiseRationalHermite:
    """Piecewise rational approximation via Hermite interpolation"""

    def __init__(self, intervals: List[Tuple[float, float]],
                 p_coeffs_list: List[np.ndarray],
                 q_coeffs_list: List[np.ndarray]):
        self.intervals = intervals
        self.p_coeffs_list = p_coeffs_list
        self.q_coeffs_list = q_coeffs_list

    @classmethod
    def from_function(cls, func: Callable, a: float, b: float, n_intervals: int,
                     deg_num: int, deg_den: int):
        """Construct piecewise rational from function"""
        K = (deg_num + deg_den - 1) // 2

        x_knots = np.linspace(a, b, n_intervals + 1)
        intervals = [(x_knots[i], x_knots[i+1]) for i in range(n_intervals)]

        p_coeffs_list = []
        q_coeffs_list = []

        for i, (ai, bi) in enumerate(intervals):
            values_a = compute_derivatives(func, ai, K)
            values_b = compute_derivatives(func, bi, K)

            p_coeffs, q_coeffs = construct_rational_hermite_robust(
                ai, bi, values_a, values_b, deg_num, deg_den
            )
            p_coeffs_list.append(p_coeffs)
            q_coeffs_list.append(q_coeffs)

        return cls(intervals, p_coeffs_list, q_coeffs_list)

    def evaluate(self, x: np.ndarray) -> np.ndarray:
        """Evaluate piecewise rational at points x"""
        result = np.zeros_like(x, dtype=float)

        for idx in range(len(x)):
            xi = x[idx]

            # Find interval
            for j, ((aj, bj), p_coeffs, q_coeffs) in enumerate(
                    zip(self.intervals, self.p_coeffs_list, self.q_coeffs_list)):

                in_interval = (aj <= xi <= bj) if j == len(self.intervals) - 1 else (aj <= xi < bj)

                if in_interval:
                    result[idx] = eval_rational(p_coeffs, q_coeffs, xi)
                    break

        return result


def test_implementation():
    """Test rational Hermite interpolation"""
    print("Testing Robust Rational Hermite Interpolation")
    print("=" * 70)

    # Test 1: Exponential
    print("\nTest 1: e^x on [0,1] with [3/2] (n=4)")
    func = lambda x: np.exp(x)
    pw = PiecewiseRationalHermite.from_function(func, 0.0, 1.0, 4, 3, 2)

    x_test = np.linspace(0, 1, 100)
    y_approx = pw.evaluate(x_test)
    y_exact = func(x_test)
    error = np.abs(y_approx - y_exact)

    print(f"  Max error: {np.nanmax(error):.2e}")
    print(f"  Mean error: {np.nanmean(error):.2e}")

    # Test 2: Sine
    print("\nTest 2: sin(x) on [0,2π] with [3/2] (n=8)")
    func = lambda x: np.sin(x)
    pw = PiecewiseRationalHermite.from_function(func, 0.0, 2*np.pi, 8, 3, 2)

    x_test = np.linspace(0, 2*np.pi, 100)
    y_approx = pw.evaluate(x_test)
    y_exact = func(x_test)
    error = np.abs(y_approx - y_exact)

    print(f"  Max error: {np.nanmax(error):.2e}")
    print(f"  Mean error: {np.nanmean(error):.2e}")

    # Test 3: Check convergence
    print("\nTest 3: Convergence study for e^x with [3/2]")
    func = lambda x: np.exp(x)
    for n in [4, 8, 16, 32]:
        pw = PiecewiseRationalHermite.from_function(func, 0.0, 1.0, n, 3, 2)
        x_test = np.linspace(0, 1, 1000)
        y_approx = pw.evaluate(x_test)
        y_exact = func(x_test)

        h = 1.0 / n
        error = y_approx - y_exact
        l2_error = np.sqrt(h * np.nansum(error**2))

        print(f"  n={n:3d}: L2 error = {l2_error:.2e}")

    print("\n" + "=" * 70)


if __name__ == "__main__":
    test_implementation()
