#!/usr/bin/env python3
"""
High-Order Finite Difference BVP Solvers

Implements 4th-order and 6th-order finite difference methods for
solving boundary value problems of the form:

    -ε*u'' + b*u' + c*u + k²*u = f(x)  on [a, b]
    u(a) = α, u(b) = β

Uses compact stencils with one-sided formulas at boundaries.
"""

import numpy as np
from scipy.linalg import solve
from dataclasses import dataclass
from typing import Callable


@dataclass
class FDResult:
    """Results from finite difference solver"""
    x: np.ndarray  # Grid points
    u: np.ndarray  # Solution values
    success: bool
    residual_norm: float = 0.0


class FourthOrderFD:
    """
    Fourth-order finite difference solver for general 2nd-order BVPs.

    Uses 5-point stencil for interior points:
        u''(x_i) ≈ [-u_{i-2} + 16u_{i-1} - 30u_i + 16u_{i+1} - u_{i+2}] / (12h²)

    One-sided stencils maintain 4th-order accuracy at boundaries.
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

    def solve(self, n: int) -> FDResult:
        """
        Solve BVP using n interior points (n+2 total points).

        Parameters
        ----------
        n : int
            Number of interior grid points

        Returns
        -------
        FDResult
            Solution at grid points
        """
        # Need at least 4 interior points for 4th-order stencil
        if n < 4:
            raise ValueError("4th-order FD requires at least n=4 interior points")

        # Create uniform grid
        x = np.linspace(self.a, self.b, n + 2)
        h = x[1] - x[0]

        # Initialize system matrix and RHS
        N = n + 2
        A = np.zeros((N, N))
        F = np.array([self.f(xi) for xi in x])

        # Build second derivative operator (D2) and first derivative operator (D1)
        D2 = np.zeros((N, N))
        D1 = np.zeros((N, N))

        # Interior points (i = 2, 3, ..., n-1) use 5-point stencil
        for i in range(2, n):
            # Second derivative: [-1, 16, -30, 16, -1] / (12h²)
            D2[i, i-2] = -1.0 / (12.0 * h**2)
            D2[i, i-1] = 16.0 / (12.0 * h**2)
            D2[i, i] = -30.0 / (12.0 * h**2)
            D2[i, i+1] = 16.0 / (12.0 * h**2)
            D2[i, i+2] = -1.0 / (12.0 * h**2)

            # First derivative: [-1, 8, 0, -8, 1] / (12h)
            D1[i, i-2] = -1.0 / (12.0 * h)
            D1[i, i-1] = 8.0 / (12.0 * h)
            D1[i, i] = 0.0
            D1[i, i+1] = -8.0 / (12.0 * h)
            D1[i, i+2] = 1.0 / (12.0 * h)

        # Near-boundary points use 4th-order one-sided stencils
        # Point i=1 (one point from left boundary)
        # Second derivative: [35, -104, 114, -56, 11] / (12h²)
        D2[1, 0] = 35.0 / (12.0 * h**2)
        D2[1, 1] = -104.0 / (12.0 * h**2)
        D2[1, 2] = 114.0 / (12.0 * h**2)
        D2[1, 3] = -56.0 / (12.0 * h**2)
        D2[1, 4] = 11.0 / (12.0 * h**2)

        # First derivative: [-11, 18, -9, 2] / (6h)
        D1[1, 0] = -11.0 / (6.0 * h)
        D1[1, 1] = 18.0 / (6.0 * h)
        D1[1, 2] = -9.0 / (6.0 * h)
        D1[1, 3] = 2.0 / (6.0 * h)

        # Point i=n (one point from right boundary)
        # Second derivative: [11, -56, 114, -104, 35] / (12h²) (reversed)
        D2[n, n-4] = 11.0 / (12.0 * h**2)
        D2[n, n-3] = -56.0 / (12.0 * h**2)
        D2[n, n-2] = 114.0 / (12.0 * h**2)
        D2[n, n-1] = -104.0 / (12.0 * h**2)
        D2[n, n] = 35.0 / (12.0 * h**2)

        # First derivative: [-2, 9, -18, 11] / (6h) (reversed)
        D1[n, n-3] = -2.0 / (6.0 * h)
        D1[n, n-2] = 9.0 / (6.0 * h)
        D1[n, n-1] = -18.0 / (6.0 * h)
        D1[n, n] = 11.0 / (6.0 * h)

        # Build operator: A = -ε*D2 + b*D1 + (c + k²)*I
        A = -self.epsilon * D2 + self.b_coeff * D1 + (self.c_coeff + self.k_squared) * np.eye(N)

        # Apply boundary conditions
        # u(x[0]) = u(a) = alpha
        A[0, :] = 0.0
        A[0, 0] = 1.0
        F[0] = self.alpha

        # u(x[n+1]) = u(b) = beta
        A[n + 1, :] = 0.0
        A[n + 1, n + 1] = 1.0
        F[n + 1] = self.beta

        # Solve linear system
        u = solve(A, F)

        # Compute residual at interior points
        residual = A @ u - F
        residual_norm = np.linalg.norm(residual[1:-1])

        return FDResult(
            x=x,
            u=u,
            success=True,
            residual_norm=residual_norm
        )


class SixthOrderFD:
    """
    Sixth-order finite difference solver for general 2nd-order BVPs.

    Uses 7-point stencil for interior points:
        u''(x_i) ≈ [2u_{i-3} - 27u_{i-2} + 270u_{i-1} - 490u_i +
                    270u_{i+1} - 27u_{i+2} + 2u_{i+3}] / (180h²)

    One-sided stencils maintain 6th-order accuracy at boundaries.
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

    def solve(self, n: int) -> FDResult:
        """
        Solve BVP using n interior points (n+2 total points).

        Parameters
        ----------
        n : int
            Number of interior grid points

        Returns
        -------
        FDResult
            Solution at grid points
        """
        # Need at least 6 interior points for 6th-order stencil
        if n < 6:
            raise ValueError("6th-order FD requires at least n=6 interior points")

        # Create uniform grid
        x = np.linspace(self.a, self.b, n + 2)
        h = x[1] - x[0]

        # Initialize system matrix and RHS
        N = n + 2
        A = np.zeros((N, N))
        F = np.array([self.f(xi) for xi in x])

        # Build second derivative operator (D2) and first derivative operator (D1)
        D2 = np.zeros((N, N))
        D1 = np.zeros((N, N))

        # Interior points (i = 3, 4, ..., n-2) use 7-point stencil
        for i in range(3, n - 1):
            # Second derivative: [2, -27, 270, -490, 270, -27, 2] / (180h²)
            D2[i, i-3] = 2.0 / (180.0 * h**2)
            D2[i, i-2] = -27.0 / (180.0 * h**2)
            D2[i, i-1] = 270.0 / (180.0 * h**2)
            D2[i, i] = -490.0 / (180.0 * h**2)
            D2[i, i+1] = 270.0 / (180.0 * h**2)
            D2[i, i+2] = -27.0 / (180.0 * h**2)
            D2[i, i+3] = 2.0 / (180.0 * h**2)

            # First derivative: [-1, 9, -45, 0, 45, -9, 1] / (60h)
            D1[i, i-3] = -1.0 / (60.0 * h)
            D1[i, i-2] = 9.0 / (60.0 * h)
            D1[i, i-1] = -45.0 / (60.0 * h)
            D1[i, i] = 0.0
            D1[i, i+1] = 45.0 / (60.0 * h)
            D1[i, i+2] = -9.0 / (60.0 * h)
            D1[i, i+3] = 1.0 / (60.0 * h)

        # Near-boundary points use 6th-order one-sided stencils
        # Point i=1 (one point from left boundary)
        # 6th-order forward difference for second derivative
        D2[1, 0] = 469.0 / (90.0 * h**2)
        D2[1, 1] = -223.0 / (10.0 * h**2)
        D2[1, 2] = 879.0 / (20.0 * h**2)
        D2[1, 3] = -949.0 / (18.0 * h**2)
        D2[1, 4] = 41.0 / (1.0 * h**2)
        D2[1, 5] = -201.0 / (10.0 * h**2)
        D2[1, 6] = 1019.0 / (180.0 * h**2)
        D2[1, 7] = -7.0 / (10.0 * h**2)

        # First derivative
        D1[1, 0] = -49.0 / (20.0 * h)
        D1[1, 1] = 6.0 / (1.0 * h)
        D1[1, 2] = -15.0 / (2.0 * h)
        D1[1, 3] = 20.0 / (3.0 * h)
        D1[1, 4] = -15.0 / (4.0 * h)
        D1[1, 5] = 6.0 / (5.0 * h)
        D1[1, 6] = -1.0 / (6.0 * h)

        # Point i=2 (two points from left boundary)
        # Mixed stencil
        D2[2, 0] = 71.0 / (90.0 * h**2)
        D2[2, 1] = -31.0 / (10.0 * h**2)
        D2[2, 2] = 31.0 / (10.0 * h**2)
        D2[2, 3] = -11.0 / (18.0 * h**2)
        D2[2, 4] = -5.0 / (1.0 * h**2)
        D2[2, 5] = 9.0 / (10.0 * h**2)
        D2[2, 6] = -11.0 / (180.0 * h**2)

        D1[2, 0] = -7.0 / (12.0 * h)
        D1[2, 1] = -1.0 / (1.0 * h)
        D1[2, 2] = 3.0 / (2.0 * h)
        D1[2, 3] = -2.0 / (3.0 * h)
        D1[2, 4] = 1.0 / (4.0 * h)
        D1[2, 5] = -1.0 / (20.0 * h)

        # Point i=n (one point from right boundary)
        # 6th-order backward difference for second derivative (mirror of i=1)
        idx = [n-7, n-6, n-5, n-4, n-3, n-2, n-1, n]
        coeffs_d2 = [-7.0/10.0, 1019.0/180.0, -201.0/10.0, 41.0/1.0,
                     -949.0/18.0, 879.0/20.0, -223.0/10.0, 469.0/90.0]
        for j, c in zip(idx, coeffs_d2):
            if 0 <= j < N:
                D2[n, j] = c / h**2

        coeffs_d1 = [1.0/6.0, -6.0/5.0, 15.0/4.0, -20.0/3.0,
                     15.0/2.0, -6.0/1.0, 49.0/20.0]
        for j, c in zip(idx[1:], coeffs_d1):
            if 0 <= j < N:
                D1[n, j] = c / h

        # Point i=n-1 (two points from right boundary)
        # Mirror of i=2
        idx = [n-6, n-5, n-4, n-3, n-2, n-1, n]
        coeffs_d2 = [-11.0/180.0, 9.0/10.0, -5.0/1.0, -11.0/18.0,
                     31.0/10.0, -31.0/10.0, 71.0/90.0]
        for j, c in zip(idx, coeffs_d2):
            if 0 <= j < N:
                D2[n - 1, j] = c / h**2

        coeffs_d1 = [1.0/20.0, -1.0/4.0, 2.0/3.0, -3.0/2.0,
                     1.0/1.0, 7.0/12.0]
        for j, c in zip(idx[1:], coeffs_d1):
            if 0 <= j < N:
                D1[n - 1, j] = c / h

        # Build operator: A = -ε*D2 + b*D1 + (c + k²)*I
        A = -self.epsilon * D2 + self.b_coeff * D1 + (self.c_coeff + self.k_squared) * np.eye(N)

        # Apply boundary conditions
        # u(x[0]) = u(a) = alpha
        A[0, :] = 0.0
        A[0, 0] = 1.0
        F[0] = self.alpha

        # u(x[n+1]) = u(b) = beta
        A[n + 1, :] = 0.0
        A[n + 1, n + 1] = 1.0
        F[n + 1] = self.beta

        # Solve linear system
        u = solve(A, F)

        # Compute residual at interior points
        residual = A @ u - F
        residual_norm = np.linalg.norm(residual[1:-1])

        return FDResult(
            x=x,
            u=u,
            success=True,
            residual_norm=residual_norm
        )
