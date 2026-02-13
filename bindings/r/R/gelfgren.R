#' Gelfgren: Piecewise Rational Interpolation
#'
#' @description
#' High-performance piecewise rational interpolation and approximation library.
#' Implements algorithms from Jan Gelfgren's 1975 paper using Bernstein polynomial
#' representation for numerical stability.
#'
#' @docType package
#' @name gelfgren-package
#' @useDynLib gelfgren, .registration = TRUE
#' @import methods
NULL

#' Create a Bernstein polynomial
#'
#' @param coeffs Numeric vector of unscaled Bernstein coefficients
#' @param a Left endpoint of interval
#' @param b Right endpoint of interval
#' @return BernsteinPolynomial object
#' @export
#' @examples
#' # Create polynomial P(x) = 1 + 2x + 3x²
#' p <- bernstein_polynomial(c(1, 2, 3), 0, 1)
#' print(p)
#'
#' # Evaluate at points
#' x <- seq(0, 1, by = 0.25)
#' y <- p$evaluate(x)
#' plot(x, y, type = "b", main = "Bernstein Polynomial")
bernstein_polynomial <- function(coeffs, a, b) {
    BernsteinPolynomial$new(coeffs, a, b)
}

#' Create a rational function
#'
#' @param numerator BernsteinPolynomial for P(x)
#' @param denominator BernsteinPolynomial for Q(x)
#' @return RationalFunction object
#' @export
#' @examples
#' num <- bernstein_polynomial(c(1, 1), 0, 1)
#' den <- bernstein_polynomial(c(1, 2), 0, 1)
#' r <- rational_function(num, den)
#'
#' x <- seq(0, 1, by = 0.1)
#' y <- r$evaluate(x)
rational_function <- function(numerator, denominator) {
    RationalFunction$new(numerator, denominator)
}

#' Create a Padé approximant
#'
#' @param coeffs Numeric vector of power series coefficients
#' @param n Numerator degree
#' @param m Denominator degree
#' @param center Expansion center
#' @param a Left endpoint of interval
#' @param b Right endpoint of interval
#' @return PadeApproximant object
#' @export
#' @examples
#' # Approximate exp(x) with [2/2] Padé
#' coeffs <- c(1, 1, 0.5, 1/6, 1/24)
#' pade <- pade_approximant(coeffs, 2, 2, 0, -1, 1)
#'
#' x <- seq(-1, 1, by = 0.25)
#' approx <- pade$evaluate(x)
#' exact <- exp(x)
#' plot(x, exact, type = "l", col = "blue", main = "Padé vs exp")
#' points(x, approx, col = "red", pch = 19)
pade_approximant <- function(coeffs, n, m, center, a, b) {
    PadeApproximant$new(coeffs, as.integer(n), as.integer(m), center, a, b)
}

#' Create a uniform mesh
#'
#' @param a Left endpoint
#' @param b Right endpoint
#' @param n Number of subintervals
#' @return Mesh object
#' @export
mesh_uniform <- function(a, b, n) {
    Mesh$uniform(a, b, as.integer(n))
}

#' Create a Chebyshev mesh
#'
#' @param a Left endpoint
#' @param b Right endpoint
#' @param n Number of subintervals
#' @return Mesh object
#' @export
mesh_chebyshev <- function(a, b, n) {
    Mesh$chebyshev(a, b, as.integer(n))
}
