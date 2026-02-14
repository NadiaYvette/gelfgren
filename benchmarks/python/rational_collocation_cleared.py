#!/usr/bin/env python3
"""
Rational Collocation BVP Solver - Cleared Form

Implements the cleared form where Q² is multiplied through before
differentiating, eliminating all divisions and providing natural pole prevention.

Theory: docs/RATIONAL_COLLOCATION_CLEARED_FORM.md
"""

import numpy as np
from scipy.optimize import least_squares
from dataclasses import dataclass
from typing import Callable, Tuple
from rational_collocation import BernsteinBasis, QConstraintType


@dataclass
class ClearedFormResult:
    """Results from cleared form rational collocation solve"""
    P_coeffs: np.ndarray      # Numerator coefficients
    Q_coeffs: np.ndarray      # Denominator coefficients
    success: bool
    message: str
    iterations: int
    residual_norm: float


class RationalCollocationCleared:
    """
    Rational collocation BVP solver using cleared form.

    Solves: -u'' = f(x) on [a,b] with u(a) = alpha, u(b) = beta
    Using rational approximant u(x) = P(x)/Q(x)

    Cleared form: Multiply ODE by Q² to eliminate divisions:
    Q² · P'' - 2Q · Q' · P' + (2Q'² - Q · Q'') · P = -Q³ · f
    """

    def __init__(self, f: Callable[[float], float],
                 a: float, b: float,
                 alpha: float, beta: float,
                 n: int, m: int,
                 num_collocation_points: int = None,
                 q_constraint: QConstraintType = QConstraintType.ENDPOINT,
                 constraint_epsilon: float = 1e-3):
        """
        Parameters
        ----------
        f : callable
            Right-hand side function f(x)
        a, b : float
            Domain endpoints
        alpha, beta : float
            Boundary values u(a) = alpha, u(b) = beta
        n, m : int
            Rational approximant degrees [n/m]
        num_collocation_points : int, optional
            Number of interior collocation points
            Default: k = n + m + 1 (for cleared form)
        q_constraint : QConstraintType, optional
            Type of constraint on Q coefficients
            Default: ENDPOINT (prevents boundary poles)
        constraint_epsilon : float, optional
            Minimum value for constrained Q coefficients
            Default: 1e-3
        """
        self.f = f
        self.a = a
        self.b = b
        self.alpha = alpha
        self.beta = beta
        self.n = n
        self.m = m
        self.q_constraint = q_constraint
        self.constraint_epsilon = constraint_epsilon

        # Determine number of collocation points
        # For cleared form, need k = n + m - 1 for square system
        # (one equation per collocation point + 2 boundary conditions)
        if num_collocation_points is None:
            self.k = n + m - 1
        else:
            self.k = num_collocation_points

        # Choose collocation points (Chebyshev nodes mapped to [a,b])
        self.collocation_points = self._chebyshev_points(self.k)

        # Total unknowns: (n+1) P coeffs + m Q coeffs
        self.num_unknowns = (n + 1) + m

        # Total equations: k collocation + 2 boundary conditions
        self.num_equations = self.k + 2

        if self.num_unknowns != self.num_equations:
            raise ValueError(
                f"System not square: {self.num_equations} equations "
                f"for {self.num_unknowns} unknowns. "
                f"Need k = n + m - 1 = {n + m - 1} collocation points."
            )

    def _chebyshev_points(self, k: int) -> np.ndarray:
        """Generate k Chebyshev points in (a,b)"""
        # Chebyshev nodes on [-1,1]
        nodes = np.cos(np.pi * (2 * np.arange(1, k+1) - 1) / (2 * k))
        # Map to [a,b]
        return self.a + (self.b - self.a) * (nodes + 1) / 2

    def _unpack_coeffs(self, coeffs: np.ndarray) -> Tuple:
        """Unpack coefficient vector into components"""
        n_p = self.n + 1
        n_q = self.m

        P_coeffs = coeffs[0:n_p]
        Q_coeffs_partial = coeffs[n_p:n_p + n_q]
        Q_coeffs = np.concatenate([[1.0], Q_coeffs_partial])  # b0 = 1 normalized

        return P_coeffs, Q_coeffs

    def _pack_coeffs(self, P_coeffs: np.ndarray, Q_coeffs: np.ndarray) -> np.ndarray:
        """Pack components into coefficient vector"""
        Q_partial = Q_coeffs[1:]  # Omit b0 = 1
        return np.concatenate([P_coeffs, Q_partial])

    def residuals(self, coeffs: np.ndarray) -> np.ndarray:
        """
        Evaluate residuals for cleared form.

        Cleared form equation at each collocation point:
        Q² · P'' - 2Q · Q' · P' + (2Q'² - Q · Q'') · P = -Q³ · f

        Plus boundary conditions:
        P(a) = Q(a) · alpha
        P(b) = Q(b) · beta
        """
        P_coeffs, Q_coeffs = self._unpack_coeffs(coeffs)

        residuals = []

        # Boundary conditions
        # At x=a: P(a) = Q(a) · alpha
        P_a = BernsteinBasis.evaluate(P_coeffs, self.a, self.a, self.b)
        Q_a = BernsteinBasis.evaluate(Q_coeffs, self.a, self.a, self.b)
        residuals.append(P_a - Q_a * self.alpha)

        # At x=b: P(b) = Q(b) · beta
        P_b = BernsteinBasis.evaluate(P_coeffs, self.b, self.a, self.b)
        Q_b = BernsteinBasis.evaluate(Q_coeffs, self.b, self.a, self.b)
        residuals.append(P_b - Q_b * self.beta)

        # Collocation equations at interior points
        for xi in self.collocation_points:
            # Evaluate P and derivatives at xi
            P = BernsteinBasis.evaluate(P_coeffs, xi, self.a, self.b)
            Px = BernsteinBasis.derivative(P_coeffs, xi, self.a, self.b, 1)
            Pxx = BernsteinBasis.derivative(P_coeffs, xi, self.a, self.b, 2)

            # Evaluate Q and derivatives at xi
            Q = BernsteinBasis.evaluate(Q_coeffs, xi, self.a, self.b)
            Qx = BernsteinBasis.derivative(Q_coeffs, xi, self.a, self.b, 1)
            Qxx = BernsteinBasis.derivative(Q_coeffs, xi, self.a, self.b, 2)

            # Cleared form: Q² · P'' - 2Q · Q' · P' + (2Q'² - Q · Q'') · P = -Q³ · f
            res = (Q**2 * Pxx -
                   2 * Q * Qx * Px +
                   (2 * Qx**2 - Q * Qxx) * P +
                   Q**3 * self.f(xi))
            residuals.append(res)

        return np.array(residuals)

    def _polynomial_initial_guess(self) -> np.ndarray:
        """Generate initial guess from polynomial finite difference solution"""
        # Solve with finite differences
        N = max(20, 2 * self.k)
        x = np.linspace(self.a, self.b, N + 1)
        h = x[1] - x[0]

        # Tridiagonal system for -u'' = f
        A = np.zeros((N + 1, N + 1))
        b = np.zeros(N + 1)

        # Interior points
        for i in range(1, N):
            A[i, i-1] = -1.0 / h**2
            A[i, i] = 2.0 / h**2
            A[i, i+1] = -1.0 / h**2
            b[i] = self.f(x[i])

        # Boundary conditions
        A[0, 0] = 1.0
        A[N, N] = 1.0
        b[0] = self.alpha
        b[N] = self.beta

        u_poly = np.linalg.solve(A, b)

        # For P, use Bernstein basis control points
        P_coeffs = np.zeros(self.n + 1)
        for i in range(self.n + 1):
            t = i / self.n if self.n > 0 else 0.5
            x_sample = self.a + t * (self.b - self.a)
            P_coeffs[i] = np.interp(x_sample, x, u_poly)

        # For Q, start with constant 1 (respecting constraints)
        if self.q_constraint in [QConstraintType.NONNEGATIVE, QConstraintType.ENDPOINT,
                                 QConstraintType.BOUNDED]:
            Q_coeffs = np.ones(self.m + 1)
        else:
            Q_coeffs = np.concatenate([[1.0], np.ones(self.m)])

        return self._pack_coeffs(P_coeffs, Q_coeffs)

    def _construct_bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        """Construct bounds on coefficients based on constraint type"""
        n_p = self.n + 1
        n_q = self.m

        # Default: no bounds
        lower = -np.inf * np.ones(self.num_unknowns)
        upper = np.inf * np.ones(self.num_unknowns)

        # Q coefficients: apply constraint based on type
        q_start = n_p
        q_end = n_p + n_q

        if self.q_constraint == QConstraintType.NONNEGATIVE:
            lower[q_start:q_end] = self.constraint_epsilon
            upper[q_start:q_end] = 10.0

        elif self.q_constraint == QConstraintType.ENDPOINT:
            lower[q_start] = self.constraint_epsilon  # b_1
            if n_q > 1:
                lower[q_end - 1] = self.constraint_epsilon  # b_m

        elif self.q_constraint == QConstraintType.BOUNDED:
            lower[q_start:q_end] = self.constraint_epsilon
            upper[q_start:q_end] = 2.0

        return lower, upper

    def solve(self, method='trf', initial_guess=None, verbose=False) -> ClearedFormResult:
        """
        Solve the BVP using cleared form rational collocation.

        Parameters
        ----------
        method : str
            Solver method: 'trf' (Trust Region Reflective) or 'lm' (Levenberg-Marquardt)
        initial_guess : ndarray, optional
            Initial guess for coefficients
        verbose : bool
            Print solver progress

        Returns
        -------
        ClearedFormResult
            Solution with coefficients and diagnostics
        """
        if initial_guess is None:
            x0 = self._polynomial_initial_guess()
        else:
            x0 = initial_guess

        # Construct bounds based on constraint type
        if self.q_constraint in [QConstraintType.NONNEGATIVE, QConstraintType.ENDPOINT,
                                 QConstraintType.BOUNDED]:
            lower_bounds, upper_bounds = self._construct_bounds()
            bounds = (lower_bounds, upper_bounds)
            solver_method = 'trf'
        else:
            bounds = (-np.inf, np.inf)
            solver_method = 'lm' if method == 'lm' else 'trf'

        # Solve
        result = least_squares(
            self.residuals,
            x0,
            method=solver_method,
            bounds=bounds,
            verbose=2 if verbose else 0,
            max_nfev=1000
        )

        # Unpack solution
        P_coeffs, Q_coeffs = self._unpack_coeffs(result.x)

        return ClearedFormResult(
            P_coeffs=P_coeffs,
            Q_coeffs=Q_coeffs,
            success=result.success,
            message=result.message,
            iterations=result.nfev,
            residual_norm=np.linalg.norm(result.fun)
        )

    def evaluate_solution(self, result: ClearedFormResult, x: np.ndarray) -> np.ndarray:
        """Evaluate rational solution u(x) = P(x)/Q(x) at points x"""
        P_vals = np.array([
            BernsteinBasis.evaluate(result.P_coeffs, xi, self.a, self.b)
            for xi in x
        ])
        Q_vals = np.array([
            BernsteinBasis.evaluate(result.Q_coeffs, xi, self.a, self.b)
            for xi in x
        ])
        return P_vals / Q_vals


if __name__ == "__main__":
    # Quick test
    print("Testing Cleared Form Rational Collocation")
    print("="*70)

    f = lambda x: 2.0
    exact = lambda x: x * (1 - x)

    solver = RationalCollocationCleared(
        f=f, a=0.0, b=1.0, alpha=0.0, beta=0.0,
        n=6, m=3,
        q_constraint=QConstraintType.ENDPOINT
    )

    print(f"Rational: [{solver.n}/{solver.m}], k={solver.k} collocation points")
    print(f"Unknowns: {solver.num_unknowns}, Equations: {solver.num_equations}")

    result = solver.solve(verbose=False)

    print(f"\nSuccess: {result.success}")
    print(f"Residual norm: {result.residual_norm:.2e}")
    print(f"Q coefficients: {result.Q_coeffs}")

    # Evaluate error
    x_test = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    u_computed = solver.evaluate_solution(result, x_test)
    u_exact = exact(x_test)
    error = np.abs(u_computed - u_exact)

    print(f"\nMax error: {np.max(error):.2e}")
    print(f"L2 error: {np.sqrt(np.mean(error**2)):.2e}")
