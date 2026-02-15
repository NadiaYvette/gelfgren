#!/usr/bin/env python3
"""
Pure Python implementation of rational Hermite interpolation

Constructs piecewise rational approximants P(x)/Q(x) by matching
function and derivatives at interval endpoints (Hermite interpolation).

Key: For [n/m] rational using k derivatives at each endpoint,
we have 2(k+1) conditions for (n+1)+(m+1) = n+m+2 unknowns.
With normalization (e.g., Q(a)=1), we have n+m+1 unknowns.
So we need: 2(k+1) = n+m+1, giving k = (n+m-1)/2

For fair comparison with Hermite polynomials:
- [3/2]: n+m=5, k=2 (match u, u', u'' at both ends) ↔ Quintic Hermite
- [5/4]: n+m=9, k=4 (match u, u', u'', u''', u'''' at both ends) ↔ Nonic Hermite
- [4/3]: n+m=7, k=3 (match u, u', u'', u''' at both ends) ↔ Septic Hermite

But we want EVEN denominator degrees for complex poles: 2, 4, 6.
So best matches:
- [3/2]: k=2 ✓
- [5/4]: k=4 ✓
- [7/6]: k=6, but need k = (7+6-1)/2 = 6 ✓
"""

import numpy as np
from scipy.optimize import least_squares, minimize
from typing import Callable, Tuple, List


def compute_derivatives(func: Callable, x: float, max_order: int) -> np.ndarray:
    """Compute function and derivatives using finite differences"""
    h = 1e-7
    derivs = [func(x)]

    if max_order >= 1:
        derivs.append((func(x + h) - func(x - h)) / (2 * h))

    if max_order >= 2:
        derivs.append((func(x + h) - 2 * func(x) + func(x - h)) / h**2)

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

    return np.array(derivs)


def rational_derivative_recursive(p_coeffs: np.ndarray, q_coeffs: np.ndarray,
                                   x: float, order: int) -> float:
    """
    Compute derivative of rational function P(x)/Q(x) recursively.

    Uses the formula: (P/Q)' = (P'Q - PQ')/Q²
    And recursion for higher orders.
    """
    def eval_poly(coeffs, x):
        return np.polyval(coeffs[::-1], x)

    def eval_poly_deriv(coeffs, x, order):
        if order == 0:
            return eval_poly(coeffs, x)
        # Derivative of polynomial
        d_coeffs = coeffs.copy()
        for _ in range(order):
            if len(d_coeffs) == 0:
                return 0.0
            d_coeffs = np.array([i * d_coeffs[i] for i in range(1, len(d_coeffs))])
        return eval_poly(d_coeffs, x)

    if order == 0:
        P = eval_poly(p_coeffs, x)
        Q = eval_poly(q_coeffs, x)
        return P / Q if abs(Q) > 1e-10 else np.nan
    elif order == 1:
        P = eval_poly(p_coeffs, x)
        Pp = eval_poly_deriv(p_coeffs, x, 1)
        Q = eval_poly(q_coeffs, x)
        Qp = eval_poly_deriv(q_coeffs, x, 1)
        return (Pp * Q - P * Qp) / Q**2 if abs(Q) > 1e-10 else np.nan
    else:
        # Use Faà di Bruno's formula or numerical differentiation
        # For simplicity, use numerical differentiation
        h = 1e-7
        f = lambda t: rational_derivative_recursive(p_coeffs, q_coeffs, t, order - 1)
        return (f(x + h) - f(x - h)) / (2 * h)


def construct_rational_hermite(a: float, b: float,
                                values_a: np.ndarray, values_b: np.ndarray,
                                deg_num: int, deg_den: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Construct rational P(x)/Q(x) matching Hermite data at endpoints.

    Parameters:
    -----------
    a, b : interval endpoints
    values_a : [u(a), u'(a), u''(a), ...]
    values_b : [u(b), u'(b), u''(b), ...]
    deg_num, deg_den : degrees of numerator and denominator

    Returns:
    --------
    p_coeffs, q_coeffs : polynomial coefficients [c0, c1, c2, ...] for sum(ci * x^i)

    Strategy: Set Q(a) = 1 (normalization), then solve nonlinear system:
        P(a)/Q(a) = u(a)
        d^k/dx^k[P(x)/Q(x)]|_{x=a} = u^(k)(a)  for k=1,...,K
        Similar at x=b
    """
    K = len(values_a) - 1  # Number of derivatives matched

    # Total unknowns: (deg_num + 1) + (deg_den + 1) - 1 = deg_num + deg_den + 1
    # (subtract 1 for normalization Q(a) = 1)
    n_unknowns = deg_num + deg_den + 1

    # Build initial guess
    # Simple strategy: P is polynomial interpolant, Q = 1
    x0 = np.zeros(n_unknowns)
    x0[:deg_num+1] = np.polyfit([a, b], [values_a[0], values_b[0]], min(1, deg_num))[::-1]
    if deg_den > 0:
        x0[deg_num+1:] = 0.0  # Q coefficients (except Q(a)=1)

    def residuals(x):
        # Extract coefficients
        p_coeffs = x[:deg_num+1]

        # Reconstruct Q coefficients with normalization Q(a) = 1
        if deg_den == 0:
            q_coeffs = np.array([1.0])
        else:
            # Q(x) = q0 + q1*x + ... = Q(a) + q1*(x-a) + q2*(x-a)^2 + ...
            # Enforce Q(a) = 1, so q0 + q1*a + q2*a^2 + ... = 1
            # We'll use shifted form: Q(x) = 1 + sum_{i=1}^m qi*(x-a)^i
            q_shifted_coeffs = np.zeros(deg_den + 1)
            q_shifted_coeffs[0] = 1.0
            q_shifted_coeffs[1:] = x[deg_num+1:]

            # Convert from (x-a) basis to x basis
            # Q(x) = 1 + q1*(x-a) + q2*(x-a)^2 + ...
            q_coeffs = np.zeros(deg_den + 1)
            for i in range(deg_den + 1):
                for j in range(i, deg_den + 1):
                    # Binomial expansion: (x-a)^j = sum_k C(j,k) x^k (-a)^(j-k)
                    from math import comb
                    q_coeffs[i] += q_shifted_coeffs[j] * comb(j, i) * (-a)**(j - i)

        res = []

        # Conditions at x = a
        for k in range(len(values_a)):
            r_deriv = rational_derivative_recursive(p_coeffs, q_coeffs, a, k)
            if np.isnan(r_deriv):
                res.append(1e10)  # Penalty for pole
            else:
                res.append(r_deriv - values_a[k])

        # Conditions at x = b
        for k in range(len(values_b)):
            r_deriv = rational_derivative_recursive(p_coeffs, q_coeffs, b, k)
            if np.isnan(r_deriv):
                res.append(1e10)
            else:
                res.append(r_deriv - values_b[k])

        # Penalty for poles in interior
        x_interior = np.linspace(a, b, 5)[1:-1]
        for xi in x_interior:
            Qi = np.polyval(q_coeffs[::-1], xi)
            if abs(Qi) < 1e-6:
                res.append(1e8 * (1e-6 - abs(Qi)))  # Penalize near-poles

        return np.array(res)

    # Solve nonlinear system
    result = least_squares(residuals, x0, method='lm', ftol=1e-10, xtol=1e-10, max_nfev=10000)

    if not result.success:
        # Try different initial guess
        x0_alt = np.random.randn(n_unknowns) * 0.1
        x0_alt[:deg_num+1] = x0[:deg_num+1]
        result = least_squares(residuals, x0_alt, method='lm', ftol=1e-10, xtol=1e-10, max_nfev=10000)

    # Extract final coefficients
    p_coeffs = result.x[:deg_num+1]

    if deg_den == 0:
        q_coeffs = np.array([1.0])
    else:
        q_shifted_coeffs = np.zeros(deg_den + 1)
        q_shifted_coeffs[0] = 1.0
        q_shifted_coeffs[1:] = result.x[deg_num+1:]

        q_coeffs = np.zeros(deg_den + 1)
        for i in range(deg_den + 1):
            for j in range(i, deg_den + 1):
                from math import comb
                q_coeffs[i] += q_shifted_coeffs[j] * comb(j, i) * (-a)**(j - i)

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
        """
        Construct piecewise rational approximation from function.

        Parameters:
        -----------
        func : function to approximate
        a, b : domain endpoints
        n_intervals : number of subintervals
        deg_num, deg_den : rational function degrees
        """
        # Determine number of derivatives to match
        K = (deg_num + deg_den - 1) // 2

        x_knots = np.linspace(a, b, n_intervals + 1)
        intervals = [(x_knots[i], x_knots[i+1]) for i in range(n_intervals)]

        p_coeffs_list = []
        q_coeffs_list = []

        for i, (ai, bi) in enumerate(intervals):
            # Compute Hermite data at endpoints
            values_a = compute_derivatives(func, ai, K)
            values_b = compute_derivatives(func, bi, K)

            # Construct rational approximant
            try:
                p_coeffs, q_coeffs = construct_rational_hermite(ai, bi, values_a, values_b,
                                                                 deg_num, deg_den)
                p_coeffs_list.append(p_coeffs)
                q_coeffs_list.append(q_coeffs)
            except Exception as e:
                print(f"Warning: Interval [{ai:.3f}, {bi:.3f}] failed: {e}")
                # Fallback to polynomial
                p_coeffs = np.polyfit([ai, bi], [values_a[0], values_b[0]], min(deg_num, 1))[::-1]
                q_coeffs = np.array([1.0])
                p_coeffs_list.append(p_coeffs)
                q_coeffs_list.append(q_coeffs)

        return cls(intervals, p_coeffs_list, q_coeffs_list)

    def evaluate(self, x: np.ndarray) -> np.ndarray:
        """Evaluate piecewise rational at points x"""
        result = np.zeros_like(x, dtype=float)

        for i in range(len(x)):
            xi = x[i]

            # Find which interval
            for j, ((aj, bj), p_coeffs, q_coeffs) in enumerate(
                    zip(self.intervals, self.p_coeffs_list, self.q_coeffs_list)):

                if j == len(self.intervals) - 1:  # Last interval
                    if aj <= xi <= bj:
                        P = np.polyval(p_coeffs[::-1], xi)
                        Q = np.polyval(q_coeffs[::-1], xi)
                        result[i] = P / Q if abs(Q) > 1e-10 else np.nan
                        break
                else:
                    if aj <= xi < bj:
                        P = np.polyval(p_coeffs[::-1], xi)
                        Q = np.polyval(q_coeffs[::-1], xi)
                        result[i] = P / Q if abs(Q) > 1e-10 else np.nan
                        break

        return result


def test_rational_hermite():
    """Test rational Hermite interpolation"""
    print("Testing Rational Hermite Interpolation")
    print("=" * 60)

    # Test 1: Exponential function with [3/2] rational
    print("\nTest 1: e^x on [0,1] with [3/2] rational (4 intervals)")
    func = lambda x: np.exp(x)

    pw_rational = PiecewiseRationalHermite.from_function(
        func, 0.0, 1.0, n_intervals=4, deg_num=3, deg_den=2
    )

    x_test = np.linspace(0, 1, 100)
    y_rational = pw_rational.evaluate(x_test)
    y_exact = func(x_test)

    error = np.abs(y_rational - y_exact)
    print(f"Max error: {np.max(error):.2e}")
    print(f"L2 error: {np.sqrt(np.mean(error**2)):.2e}")

    # Test 2: Sine function
    print("\nTest 2: sin(x) on [0,2π] with [3/2] rational (8 intervals)")
    func = lambda x: np.sin(x)

    pw_rational = PiecewiseRationalHermite.from_function(
        func, 0.0, 2*np.pi, n_intervals=8, deg_num=3, deg_den=2
    )

    x_test = np.linspace(0, 2*np.pi, 100)
    y_rational = pw_rational.evaluate(x_test)
    y_exact = func(x_test)

    error = np.abs(y_rational - y_exact)
    print(f"Max error: {np.max(error):.2e}")
    print(f"L2 error: {np.sqrt(np.mean(error**2)):.2e}")

    print("\n" + "=" * 60)
    print("Tests complete!")


if __name__ == "__main__":
    test_rational_hermite()
