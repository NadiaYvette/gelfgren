#!/usr/bin/env python3
"""
Spectral Collocation BVP Solvers

Implements Chebyshev and Legendre spectral collocation methods for
solving boundary value problems of the form:

    -ε*u'' + b*u' + c*u + k²*u = f(x)  on [a, b]
    u(a) = α, u(b) = β

Uses differentiation matrices at Gauss-Lobatto collocation points.
"""

import numpy as np
from scipy.linalg import solve
from dataclasses import dataclass
from typing import Callable, Optional


@dataclass
class SpectralResult:
    """Results from spectral collocation solver"""
    x: np.ndarray  # Collocation points (physical domain)
    u: np.ndarray  # Solution values
    success: bool
    residual_norm: float = 0.0


class ChebyshevSpectralSolver:
    """
    Chebyshev spectral collocation solver for general 2nd-order BVPs.

    Uses Chebyshev-Gauss-Lobatto points and differentiation matrix
    via the Trefethen algorithm (Spectral Methods in MATLAB, 2000).
    """

    def __init__(self, f: Callable, a: float, b: float,
                 alpha: float, beta: float,
                 epsilon: float = 1.0, b_coeff: float = 0.0,
                 c_coeff: float = 0.0, k_squared: float = 0.0):
        """
        Parameters
        ----------
        f : callable
            Forcing function f(x)
        a, b : float
            Domain boundaries
        alpha, beta : float
            Boundary values u(a) = alpha, u(b) = beta
        epsilon : float
            Diffusion coefficient (default 1.0)
        b_coeff : float
            Advection coefficient (default 0.0)
        c_coeff : float
            Reaction coefficient (default 0.0)
        k_squared : float
            Helmholtz parameter (default 0.0)
        """
        self.f = f
        self.a = a
        self.b = b
        self.alpha = alpha
        self.beta = beta
        self.epsilon = epsilon
        self.b_coeff = b_coeff
        self.c_coeff = c_coeff
        self.k_squared = k_squared

    def chebyshev_points(self, N: int) -> np.ndarray:
        """
        Generate N+1 Chebyshev-Gauss-Lobatto points in [-1, 1].

        Points are: x_j = cos(π*j/N) for j = 0, 1, ..., N
        """
        return np.cos(np.pi * np.arange(N + 1) / N)

    def chebyshev_diff_matrix(self, N: int) -> np.ndarray:
        """
        Chebyshev differentiation matrix via Trefethen algorithm.

        Computes first derivative matrix D such that:
            u'(x_j) ≈ Σ_k D_jk * u(x_k)

        Reference: Trefethen, "Spectral Methods in MATLAB", SIAM 2000
        """
        if N == 0:
            return np.zeros((1, 1))

        x = self.chebyshev_points(N)

        # Weights for each point
        c = np.ones(N + 1)
        c[0] = 2.0
        c[N] = 2.0
        c = c * (-1.0) ** np.arange(N + 1)

        # Build differentiation matrix
        X = np.tile(x.reshape(-1, 1), (1, N + 1))
        dX = X - X.T

        D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))  # off-diagonal entries
        D = D - np.diag(np.sum(D, axis=1))  # diagonal entries

        return D

    def transform_to_physical(self, x_ref: np.ndarray) -> np.ndarray:
        """Transform from reference domain [-1, 1] to physical domain [a, b]."""
        return 0.5 * (self.b - self.a) * (x_ref + 1) + self.a

    def solve(self, N: int) -> SpectralResult:
        """
        Solve BVP using N+1 Chebyshev-Gauss-Lobatto points.

        Parameters
        ----------
        N : int
            Polynomial degree (N+1 collocation points)

        Returns
        -------
        SpectralResult
            Solution at collocation points
        """
        # Get collocation points in reference domain [-1, 1]
        x_ref = self.chebyshev_points(N)

        # Transform to physical domain [a, b]
        # Note: Chebyshev points are in decreasing order (1 → -1)
        # We reverse them to get increasing order (a → b)
        x_ref = x_ref[::-1]
        x = self.transform_to_physical(x_ref)

        # Get differentiation matrices in reference domain
        D = self.chebyshev_diff_matrix(N)[::-1, ::-1]  # reverse for increasing order

        # Scale for physical domain: d/dx = (2/(b-a)) * d/dξ
        scale = 2.0 / (self.b - self.a)
        D1 = scale * D
        D2 = scale**2 * (D @ D)

        # Build linear system: L*u = F
        # Operator: L = -ε*D² + b*D + (c + k²)*I
        L = -self.epsilon * D2 + self.b_coeff * D1 + (self.c_coeff + self.k_squared) * np.eye(N + 1)

        # Right-hand side
        F = np.array([self.f(xi) for xi in x])

        # Apply boundary conditions
        # u(x[0]) = u(a) = alpha
        L[0, :] = 0.0
        L[0, 0] = 1.0
        F[0] = self.alpha

        # u(x[N]) = u(b) = beta
        L[N, :] = 0.0
        L[N, N] = 1.0
        F[N] = self.beta

        # Solve linear system
        u = solve(L, F)

        # Compute residual at interior points
        residual = L @ u - F
        residual_norm = np.linalg.norm(residual[1:-1])

        return SpectralResult(
            x=x,
            u=u,
            success=True,
            residual_norm=residual_norm
        )


class LegendreSpectralSolver:
    """
    Legendre spectral collocation solver for general 2nd-order BVPs.

    Uses Legendre-Gauss-Lobatto points and differentiation matrix.
    """

    def __init__(self, f: Callable, a: float, b: float,
                 alpha: float, beta: float,
                 epsilon: float = 1.0, b_coeff: float = 0.0,
                 c_coeff: float = 0.0, k_squared: float = 0.0):
        """
        Parameters
        ----------
        f : callable
            Forcing function f(x)
        a, b : float
            Domain boundaries
        alpha, beta : float
            Boundary values u(a) = alpha, u(b) = beta
        epsilon : float
            Diffusion coefficient (default 1.0)
        b_coeff : float
            Advection coefficient (default 0.0)
        c_coeff : float
            Reaction coefficient (default 0.0)
        k_squared : float
            Helmholtz parameter (default 0.0)
        """
        self.f = f
        self.a = a
        self.b = b
        self.alpha = alpha
        self.beta = beta
        self.epsilon = epsilon
        self.b_coeff = b_coeff
        self.c_coeff = c_coeff
        self.k_squared = k_squared

    def legendre_gauss_lobatto_points(self, N: int) -> np.ndarray:
        """
        Compute N+1 Legendre-Gauss-Lobatto points in [-1, 1].

        These are zeros of (1-x²)*P'_N(x) where P_N is Legendre polynomial.
        Uses Newton iteration.
        """
        if N == 0:
            return np.array([0.0])
        if N == 1:
            return np.array([-1.0, 1.0])

        # Initial guess using Chebyshev points
        x = np.cos(np.pi * np.arange(N + 1) / N)

        # Newton iteration to find LGL points
        P = np.zeros((N + 1, N + 2))

        for _ in range(10):  # Usually converges in 3-4 iterations
            # Compute Legendre polynomial P_N and its derivative at x
            P[:, 0] = 1.0
            P[:, 1] = x

            for k in range(2, N + 2):
                P[:, k] = ((2*k - 1) * x * P[:, k-1] - (k - 1) * P[:, k-2]) / k

            # Function: f(x) = (1-x²)*P'_N(x)
            # Derivative: f'(x) = -2x*P'_N(x) + (1-x²)*P''_N(x)
            # P'_N(x) = N/(1-x²) * (P_{N-1}(x) - x*P_N(x))

            x_old = x.copy()

            # Update (avoid division by zero at endpoints)
            for j in range(1, N):
                P_N = P[j, N]
                P_N_minus_1 = P[j, N - 1]

                # Derivative of P_N
                dP_N = N / (1 - x[j]**2) * (P_N_minus_1 - x[j] * P_N)

                # Second derivative of P_N
                ddP_N = (2 * x[j] * dP_N - N * (N + 1) * P_N) / (1 - x[j]**2)

                # Newton update for (1-x²)*P'_N(x) = 0
                f = (1 - x[j]**2) * dP_N
                df = -2 * x[j] * dP_N + (1 - x[j]**2) * ddP_N

                x[j] = x[j] - f / df

            # Endpoints are fixed
            x[0] = 1.0
            x[N] = -1.0

            if np.max(np.abs(x - x_old)) < 1e-15:
                break

        return x[::-1]  # Return in increasing order

    def legendre_diff_matrix(self, N: int, x: np.ndarray) -> np.ndarray:
        """
        Legendre differentiation matrix at points x.

        Uses barycentric formula for differentiation matrix.
        """
        if N == 0:
            return np.zeros((1, 1))

        # Compute Legendre polynomial P_N at points x
        P = np.zeros((N + 1, N + 2))
        P[:, 0] = 1.0
        P[:, 1] = x

        for k in range(2, N + 2):
            P[:, k] = ((2*k - 1) * x * P[:, k-1] - (k - 1) * P[:, k-2]) / k

        P_N = P[:, N]

        # Build differentiation matrix
        D = np.zeros((N + 1, N + 1))

        for i in range(N + 1):
            for j in range(N + 1):
                if i != j:
                    D[i, j] = (P_N[i] / P_N[j]) / (x[i] - x[j])
                else:
                    if i == 0 or i == N:
                        # Endpoints
                        D[i, i] = x[i] / (1 - x[i]**2)
                    else:
                        # Interior points: diagonal is zero for LGL
                        D[i, i] = 0.0

        # Correct diagonal at endpoints
        D[0, 0] = -N * (N + 1) / 4.0
        D[N, N] = N * (N + 1) / 4.0

        return D

    def transform_to_physical(self, x_ref: np.ndarray) -> np.ndarray:
        """Transform from reference domain [-1, 1] to physical domain [a, b]."""
        return 0.5 * (self.b - self.a) * (x_ref + 1) + self.a

    def solve(self, N: int) -> SpectralResult:
        """
        Solve BVP using N+1 Legendre-Gauss-Lobatto points.

        Parameters
        ----------
        N : int
            Polynomial degree (N+1 collocation points)

        Returns
        -------
        SpectralResult
            Solution at collocation points
        """
        # Get collocation points in reference domain [-1, 1]
        x_ref = self.legendre_gauss_lobatto_points(N)

        # Transform to physical domain [a, b]
        x = self.transform_to_physical(x_ref)

        # Get differentiation matrix in reference domain
        D = self.legendre_diff_matrix(N, x_ref)

        # Scale for physical domain: d/dx = (2/(b-a)) * d/dξ
        scale = 2.0 / (self.b - self.a)
        D1 = scale * D
        D2 = scale**2 * (D @ D)

        # Build linear system: L*u = F
        # Operator: L = -ε*D² + b*D + (c + k²)*I
        L = -self.epsilon * D2 + self.b_coeff * D1 + (self.c_coeff + self.k_squared) * np.eye(N + 1)

        # Right-hand side
        F = np.array([self.f(xi) for xi in x])

        # Apply boundary conditions
        # u(x[0]) = u(a) = alpha
        L[0, :] = 0.0
        L[0, 0] = 1.0
        F[0] = self.alpha

        # u(x[N]) = u(b) = beta
        L[N, :] = 0.0
        L[N, N] = 1.0
        F[N] = self.beta

        # Solve linear system
        u = solve(L, F)

        # Compute residual at interior points
        residual = L @ u - F
        residual_norm = np.linalg.norm(residual[1:-1])

        return SpectralResult(
            x=x,
            u=u,
            success=True,
            residual_norm=residual_norm
        )
