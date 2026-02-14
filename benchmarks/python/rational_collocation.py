#!/usr/bin/env python3
"""
Rational Collocation BVP Solver - Quadratic Formulation

Implements the quadratic formulation where u(x_i) and u'(x_i) are
explicit unknowns, resulting in a bilinear (quadratic) system.

Theory: docs/RATIONAL_COLLOCATION_QUADRATIC_FORM.md
"""

import numpy as np
from scipy.optimize import fsolve, least_squares
from dataclasses import dataclass
from typing import Callable, Tuple, List
from enum import Enum


class QConstraintType(Enum):
    """Types of constraints on Q coefficients to prevent poles"""
    NONE = "none"  # No constraints (may find spurious poles)
    REGULARIZATION = "regularization"  # Penalty term pushing Q toward constant
    NONNEGATIVE = "nonnegative"  # All Q coefficients >= 0 (prevents all poles)
    ENDPOINT = "endpoint"  # Only b_0, b_m >= epsilon (prevents boundary poles)
    BOUNDED = "bounded"  # Q coefficients in [0, 2] (non-negative + bounded)


@dataclass
class CollocationResult:
    """Results from rational collocation solve"""
    P_coeffs: np.ndarray      # Numerator coefficients
    Q_coeffs: np.ndarray      # Denominator coefficients
    u_values: np.ndarray      # Solution values at collocation points
    up_values: np.ndarray     # Solution derivatives at collocation points
    success: bool
    message: str
    iterations: int
    residual_norm: float


class BernsteinBasis:
    """Helper for Bernstein polynomial evaluation"""

    @staticmethod
    def binomial(n: int, k: int) -> float:
        """Binomial coefficient C(n,k)"""
        if k > n or k < 0:
            return 0.0
        if k == 0 or k == n:
            return 1.0
        result = 1.0
        for i in range(min(k, n - k)):
            result *= (n - i) / (i + 1)
        return result

    @staticmethod
    def evaluate(coeffs: np.ndarray, x: float, a: float, b: float) -> float:
        """Evaluate Bernstein polynomial at x ∈ [a,b]"""
        n = len(coeffs) - 1
        t = (x - a) / (b - a)  # Normalize to [0,1]

        result = 0.0
        for i, c in enumerate(coeffs):
            B_i = BernsteinBasis.binomial(n, i) * t**i * (1-t)**(n-i)
            result += c * B_i

        return result

    @staticmethod
    def derivative(coeffs: np.ndarray, x: float, a: float, b: float, order: int = 1) -> float:
        """Evaluate derivative of Bernstein polynomial"""
        if order == 0:
            return BernsteinBasis.evaluate(coeffs, x, a, b)

        n = len(coeffs) - 1
        delta = b - a

        # First derivative coefficients
        if order >= 1:
            if n == 0:
                deriv_coeffs = np.array([0.0])
            else:
                deriv_coeffs = n / delta * np.diff(coeffs)

        # Higher derivatives
        for _ in range(1, order):
            if len(deriv_coeffs) == 0:
                return 0.0
            n_curr = len(deriv_coeffs) - 1
            if n_curr == 0:
                deriv_coeffs = np.array([0.0])
            else:
                deriv_coeffs = n_curr / delta * np.diff(deriv_coeffs)

        if len(deriv_coeffs) == 0 or (len(deriv_coeffs) == 1 and order > 0):
            return 0.0

        return BernsteinBasis.evaluate(deriv_coeffs, x, a, b)


class RationalCollocationQuadratic:
    """
    Rational collocation BVP solver using quadratic formulation.

    Solves: -u'' = f(x) on [a,b] with u(a) = alpha, u(b) = beta
    Using rational approximant u(x) = P(x)/Q(x)

    Formulation: Treat u(x_i) and u'(x_i) as explicit unknowns,
    resulting in quadratic nonlinearity with bilinear structure.
    """

    def __init__(self, f: Callable[[float], float],
                 a: float, b: float,
                 alpha: float, beta: float,
                 n: int, m: int,
                 num_collocation_points: int = None,
                 q_constraint: QConstraintType = QConstraintType.NONNEGATIVE,
                 q_regularization: float = 0.0,
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
            Default: k = n + m - 1 (for well-determined system)
        q_constraint : QConstraintType, optional
            Type of constraint on Q coefficients to prevent poles
            Default: NONNEGATIVE (all Q coeffs >= 0)
        q_regularization : float, optional
            Regularization weight (only used if q_constraint=REGULARIZATION)
            Default: 0.0
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
        self.q_regularization = q_regularization
        self.constraint_epsilon = constraint_epsilon

        # Determine number of collocation points
        if num_collocation_points is None:
            self.k = n + m - 1
        else:
            self.k = num_collocation_points

        # Choose collocation points (Chebyshev nodes mapped to [a,b])
        self.collocation_points = self._chebyshev_points(self.k)

        # Total unknowns: (n+1) P coeffs + m Q coeffs + 2k (u and u' values)
        self.num_unknowns = (n + 1) + m + 2 * self.k

        # Total equations: 3k collocation + 2 boundary conditions
        self.num_equations = 3 * self.k + 2

        # Add regularization terms if using that constraint type
        if self.q_constraint == QConstraintType.REGULARIZATION and self.q_regularization > 0:
            self.num_equations += m  # One regularization term per Q coefficient

        # For overdetermined systems (with regularization), use least squares
        # For square systems (with inequality constraints), check dimensions
        using_regularization = (self.q_constraint == QConstraintType.REGULARIZATION and
                               self.q_regularization > 0)
        if not using_regularization and self.num_unknowns != self.num_equations:
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

        u_and_up = coeffs[n_p + n_q:]
        u_vals = u_and_up[0::2]  # Even indices
        up_vals = u_and_up[1::2]  # Odd indices

        return P_coeffs, Q_coeffs, u_vals, up_vals

    def _pack_coeffs(self, P_coeffs: np.ndarray, Q_coeffs: np.ndarray,
                     u_vals: np.ndarray, up_vals: np.ndarray) -> np.ndarray:
        """Pack components into coefficient vector"""
        Q_partial = Q_coeffs[1:]  # Omit b0 = 1
        u_and_up = np.empty(2 * len(u_vals))
        u_and_up[0::2] = u_vals
        u_and_up[1::2] = up_vals
        return np.concatenate([P_coeffs, Q_partial, u_and_up])

    def residuals(self, coeffs: np.ndarray) -> np.ndarray:
        """
        Evaluate residuals for quadratic formulation.

        Equations:
        (1) P(x_i) - Q(x_i)·u(x_i) = 0
        (2) P'(x_i) - Q'(x_i)·u(x_i) - Q(x_i)·u'(x_i) = 0
        (3) P''(x_i) - Q''(x_i)·u(x_i) - 2Q'(x_i)·u'(x_i) + Q(x_i)·f(x_i) = 0

        Plus boundary conditions:
        u(a) = alpha (enforced on u_values)
        u(b) = beta  (enforced on u_values)
        """
        P_coeffs, Q_coeffs, u_vals, up_vals = self._unpack_coeffs(coeffs)

        residuals = []

        # Boundary conditions
        # At x=a: P(a) = Q(a)·alpha (enforce exactly)
        P_a = BernsteinBasis.evaluate(P_coeffs, self.a, self.a, self.b)
        Q_a = BernsteinBasis.evaluate(Q_coeffs, self.a, self.a, self.b)
        residuals.append(P_a - Q_a * self.alpha)

        # At x=b: P(b) = Q(b)·beta
        P_b = BernsteinBasis.evaluate(P_coeffs, self.b, self.a, self.b)
        Q_b = BernsteinBasis.evaluate(Q_coeffs, self.b, self.a, self.b)
        residuals.append(P_b - Q_b * self.beta)

        # Collocation equations at interior points
        for i, xi in enumerate(self.collocation_points):
            ui = u_vals[i]
            upi = up_vals[i]

            # Evaluate P, P', P'' at xi
            P = BernsteinBasis.evaluate(P_coeffs, xi, self.a, self.b)
            Px = BernsteinBasis.derivative(P_coeffs, xi, self.a, self.b, 1)
            Pxx = BernsteinBasis.derivative(P_coeffs, xi, self.a, self.b, 2)

            # Evaluate Q, Q', Q'' at xi
            Q = BernsteinBasis.evaluate(Q_coeffs, xi, self.a, self.b)
            Qx = BernsteinBasis.derivative(Q_coeffs, xi, self.a, self.b, 1)
            Qxx = BernsteinBasis.derivative(Q_coeffs, xi, self.a, self.b, 2)

            # Equation (1): P(xi) = Q(xi)·u(xi)
            res1 = P - Q * ui
            residuals.append(res1)

            # Equation (2): P'(xi) = Q'(xi)·u(xi) + Q(xi)·u'(xi)
            res2 = Px - Qx * ui - Q * upi
            residuals.append(res2)

            # Equation (3): P''(xi) = Q''(xi)·u(xi) + 2Q'(xi)·u'(xi) - Q(xi)·f(xi)
            res3 = Pxx - Qxx * ui - 2 * Qx * upi + Q * self.f(xi)
            residuals.append(res3)

        # Add regularization terms if using that constraint type
        # For Bernstein polynomials, constant Q=1 means all coeffs = 1
        if self.q_constraint == QConstraintType.REGULARIZATION and self.q_regularization > 0:
            Q_partial = Q_coeffs[1:]  # Coefficients other than b0=1
            for q_coeff in Q_partial:
                # Penalize deviation from 1 (not from 0!)
                residuals.append(self.q_regularization * (q_coeff - 1.0))

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

        # Interpolate to collocation points
        u_vals = np.interp(self.collocation_points, x, u_poly)

        # Estimate derivatives
        up_vals = np.gradient(u_poly, h)
        up_vals = np.interp(self.collocation_points, x, up_vals)

        # Initial guess for P and Q
        # For P, use Bernstein basis control points
        # Sample u_poly at Bernstein control point locations
        P_coeffs = np.zeros(self.n + 1)
        for i in range(self.n + 1):
            # Bernstein control point i is at position i/n in parameter space
            t = i / self.n if self.n > 0 else 0.5
            x_sample = self.a + t * (self.b - self.a)
            P_coeffs[i] = np.interp(x_sample, x, u_poly)

        # For Q, initial guess depends on constraint type
        if self.q_constraint in [QConstraintType.NONNEGATIVE, QConstraintType.ENDPOINT,
                                 QConstraintType.BOUNDED]:
            # Need all coefficients to be positive and close to 1 for constant
            # Start with Q ≈ 1 everywhere, respecting the constraint
            Q_coeffs = np.ones(self.m + 1)
        else:
            # No constraints or regularization: start with Q = constant 1
            Q_coeffs = np.concatenate([[1.0], np.ones(self.m)])

        return self._pack_coeffs(P_coeffs, Q_coeffs, u_vals, up_vals)

    def _construct_bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Construct bounds on coefficients based on constraint type.

        Returns
        -------
        lower_bounds, upper_bounds : ndarray
            Lower and upper bounds for each coefficient
        """
        n_p = self.n + 1  # Number of P coefficients
        n_q = self.m      # Number of Q coefficients (excluding b0=1)
        n_u = 2 * self.k  # Number of u and u' values

        # Default: no bounds
        lower = -np.inf * np.ones(self.num_unknowns)
        upper = np.inf * np.ones(self.num_unknowns)

        # P coefficients: no constraints
        # (Could add bounds here if needed for stability)

        # Q coefficients: apply constraint based on type
        q_start = n_p
        q_end = n_p + n_q

        if self.q_constraint == QConstraintType.NONNEGATIVE:
            # All Q coefficients >= epsilon
            lower[q_start:q_end] = self.constraint_epsilon
            # Also allow reasonable upper bound to prevent extreme values
            upper[q_start:q_end] = 10.0

        elif self.q_constraint == QConstraintType.ENDPOINT:
            # Only constrain first and last Q coefficients
            # Note: Q_coeffs_partial[0] corresponds to b_1, Q_coeffs_partial[m-1] to b_m
            # (b_0 = 1 is fixed and not in the unknowns)
            lower[q_start] = self.constraint_epsilon  # b_1
            if n_q > 1:
                lower[q_end - 1] = self.constraint_epsilon  # b_m

        elif self.q_constraint == QConstraintType.BOUNDED:
            # Q coefficients in [epsilon, 2]
            lower[q_start:q_end] = self.constraint_epsilon
            upper[q_start:q_end] = 2.0

        # u and u' values: very permissive bounds
        # Don't over-constrain the solution values
        u_start = q_end
        lower[u_start:] = -1000.0  # Very permissive
        upper[u_start:] = 1000.0

        return lower, upper

    def solve(self, method='lm', initial_guess=None, verbose=False) -> CollocationResult:
        """
        Solve the BVP using rational collocation.

        Parameters
        ----------
        method : str
            Solver method: 'lm' (Levenberg-Marquardt), 'trf' (Trust Region Reflective)
        initial_guess : ndarray, optional
            Initial guess for coefficients
        verbose : bool
            Print solver progress

        Returns
        -------
        CollocationResult
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
            # Use 'trf' method which supports bounds
            solver_method = 'trf'
        else:
            bounds = (-np.inf, np.inf)
            # Use 'lm' for unconstrained (faster)
            solver_method = 'lm' if method == 'lm' else 'trf'

        if method in ['lm', 'trf']:
            # Levenberg-Marquardt or Trust Region Reflective
            result = least_squares(
                self.residuals,
                x0,
                method=solver_method,
                bounds=bounds,
                verbose=2 if verbose else 0,
                max_nfev=1000  # Increase max iterations for constrained problems
            )
            success = result.success
            message = result.message
            iterations = result.nfev
            residual_norm = np.linalg.norm(result.fun)
            solution = result.x

        elif method == 'fsolve':
            # Newton's method (no bounds support)
            if self.q_constraint != QConstraintType.NONE and self.q_constraint != QConstraintType.REGULARIZATION:
                raise ValueError("fsolve does not support bounds. Use method='trf' with constraints.")

            solution, info, ier, msg = fsolve(
                self.residuals,
                x0,
                full_output=True
            )
            success = (ier == 1)
            message = msg
            iterations = info['nfev']
            residual_norm = np.linalg.norm(info['fvec'])

        else:
            raise ValueError(f"Unknown method: {method}. Use 'lm', 'trf', or 'fsolve'.")

        # Unpack solution
        P_coeffs, Q_coeffs, u_vals, up_vals = self._unpack_coeffs(solution)

        return CollocationResult(
            P_coeffs=P_coeffs,
            Q_coeffs=Q_coeffs,
            u_values=u_vals,
            up_values=up_vals,
            success=success,
            message=message,
            iterations=iterations,
            residual_norm=residual_norm
        )

    def evaluate_solution(self, result: CollocationResult, x: np.ndarray) -> np.ndarray:
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


def test_simple_poisson():
    """Test on simple problem: -u'' = 2, u(0) = u(1) = 0"""
    print("="*70)
    print("Test: -u'' = 2 on [0,1], u(0) = u(1) = 0")
    print("Exact: u(x) = x(1-x)")
    print("="*70)

    # Problem setup
    f = lambda x: 2.0
    exact = lambda x: x * (1 - x)

    # Test different constraint types
    constraint_types = [
        (QConstraintType.NONNEGATIVE, "Non-negative Q coefficients"),
        (QConstraintType.ENDPOINT, "Endpoint constraints only"),
        (QConstraintType.BOUNDED, "Bounded Q ∈ [ε, 2]"),
        (QConstraintType.REGULARIZATION, "Regularization (penalty)"),
    ]

    results = []

    for constraint_type, description in constraint_types:
        print(f"\n{'-'*70}")
        print(f"Constraint: {description}")
        print(f"{'-'*70}")

        # Solve with [6/3] rational
        solver = RationalCollocationQuadratic(
            f=f, a=0.0, b=1.0, alpha=0.0, beta=0.0,
            n=6, m=3,
            q_constraint=constraint_type,
            q_regularization=100.0 if constraint_type == QConstraintType.REGULARIZATION else 0.0,
            constraint_epsilon=1e-3
        )

        print(f"Rational: [{solver.n}/{solver.m}], k={solver.k} collocation points")

        try:
            # Solve
            result = solver.solve(method='trf', verbose=False)

            print(f"Success: {result.success}")
            print(f"Iterations: {result.iterations}")
            print(f"Residual norm: {result.residual_norm:.2e}")
            print(f"Q coefficients: {result.Q_coeffs}")

            # Check Q values
            x_check = np.linspace(0, 1, 100)
            Q_vals = np.array([BernsteinBasis.evaluate(result.Q_coeffs, xi, 0.0, 1.0)
                              for xi in x_check])
            print(f"Q range: [{np.min(Q_vals):.6f}, {np.max(Q_vals):.6f}]")

            # Evaluate error
            x_test = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
            u_computed = solver.evaluate_solution(result, x_test)
            u_exact = exact(x_test)
            error = np.abs(u_computed - u_exact)

            print(f"Max error: {np.max(error):.2e}")
            print(f"L2 error: {np.sqrt(np.mean(error**2)):.2e}")

            results.append((constraint_type, result, np.max(error)))

        except Exception as e:
            print(f"FAILED: {str(e)[:100]}")
            results.append((constraint_type, None, float('inf')))

    print(f"\n{'='*70}")
    print("Summary")
    print(f"{'='*70}")
    print(f"{'Constraint':<35} {'Max Error':<15} {'Q range':<15}")
    print(f"{'-'*70}")
    for constraint_type, result, max_error in results:
        if result is not None:
            x_check = np.linspace(0, 1, 100)
            Q_vals = np.array([BernsteinBasis.evaluate(result.Q_coeffs, xi, 0.0, 1.0)
                              for xi in x_check])
            q_range = f"[{np.min(Q_vals):.3f}, {np.max(Q_vals):.3f}]"
        else:
            q_range = "FAILED"

        print(f"{constraint_type.value:<35} {max_error:<15.2e} {q_range:<15}")

    return results


if __name__ == "__main__":
    test_simple_poisson()
