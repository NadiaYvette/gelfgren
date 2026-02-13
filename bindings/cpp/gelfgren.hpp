/**
 * @file gelfgren.hpp
 * @brief C++ RAII wrappers for Gelfgren C API
 *
 * Provides safe, object-oriented interface to Gelfgren library with automatic
 * memory management, exception-based error handling, and operator overloading.
 */

#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

extern "C" {
#include "../../include/gelfgren.h"
}

namespace gelfgren {

/**
 * @brief Exception thrown when Gelfgren operations fail
 */
class GelfgrenException : public std::runtime_error {
public:
    explicit GelfgrenException(const std::string& message)
        : std::runtime_error(message) {}

    explicit GelfgrenException(GelfgrenErrorCode code)
        : std::runtime_error(get_error_message(code)) {}

private:
    static std::string get_error_message(GelfgrenErrorCode code) {
        const char* msg = gelfgren_last_error_message();
        if (msg != nullptr) {
            return std::string(msg);
        }
        return "Unknown error (code: " + std::to_string(static_cast<int>(code)) + ")";
    }
};

/**
 * @brief Check error code and throw exception if not successful
 */
inline void check_error(GelfgrenErrorCode code) {
    if (code != GelfgrenErrorCode::SUCCESS) {
        throw GelfgrenException(code);
    }
}

/**
 * @brief Check pointer and throw exception if null
 */
template<typename T>
inline T* check_ptr(T* ptr) {
    if (ptr == nullptr) {
        throw GelfgrenException("Null pointer returned from Gelfgren function");
    }
    return ptr;
}

// Forward declarations
class RationalFunction;

/**
 * @brief RAII wrapper for Bernstein polynomial
 *
 * Provides automatic memory management and C++ operator overloading.
 * Bernstein polynomials use the numerically stable Bernstein basis.
 *
 * @example
 * ```cpp
 * // Create polynomial P(x) = 1 + 2x + 3x²
 * BernsteinPolynomial p({1.0, 2.0, 3.0}, 0.0, 1.0);
 * double value = p(0.5);  // Evaluate at x = 0.5
 * auto derivative = p.derivative();
 * ```
 */
class BernsteinPolynomial {
public:
    /**
     * @brief Construct from unscaled coefficients
     * @param coeffs Unscaled Bernstein coefficients [C₀, C₁, ..., Cₙ]
     * @param a Left endpoint of interval
     * @param b Right endpoint of interval
     */
    BernsteinPolynomial(const std::vector<double>& coeffs, double a, double b)
        : ptr_(check_ptr(gelfgren_bernstein_create(
            coeffs.data(), coeffs.size() - 1, a, b)),
            &gelfgren_bernstein_free) {}

    /**
     * @brief Evaluate polynomial at point x
     * @param x Evaluation point
     * @return P(x)
     */
    double operator()(double x) const {
        double result;
        check_error(gelfgren_bernstein_evaluate(ptr_.get(), x, &result));
        return result;
    }

    /**
     * @brief Evaluate polynomial at point x (named method)
     */
    double evaluate(double x) const {
        return (*this)(x);
    }

    /**
     * @brief Get polynomial degree
     */
    size_t degree() const {
        size_t result;
        check_error(gelfgren_bernstein_degree(ptr_.get(), &result));
        return result;
    }

    /**
     * @brief Get interval [a, b]
     */
    std::pair<double, double> interval() const {
        double a, b;
        check_error(gelfgren_bernstein_interval(ptr_.get(), &a, &b));
        return {a, b};
    }

    /**
     * @brief Compute derivative polynomial
     */
    BernsteinPolynomial derivative() const {
        return BernsteinPolynomial(
            check_ptr(gelfgren_bernstein_derivative(ptr_.get())));
    }

    /**
     * @brief Compute integral (antiderivative) polynomial
     */
    BernsteinPolynomial integral() const {
        return BernsteinPolynomial(
            check_ptr(gelfgren_bernstein_integral(ptr_.get())));
    }

    /**
     * @brief Elevate polynomial degree by 1
     */
    BernsteinPolynomial elevate() const {
        return BernsteinPolynomial(
            check_ptr(gelfgren_bernstein_elevate(ptr_.get())));
    }

    /**
     * @brief Add two polynomials
     */
    BernsteinPolynomial operator+(const BernsteinPolynomial& other) const {
        return BernsteinPolynomial(
            check_ptr(gelfgren_bernstein_add(ptr_.get(), other.ptr_.get())));
    }

    /**
     * @brief Subtract two polynomials
     */
    BernsteinPolynomial operator-(const BernsteinPolynomial& other) const {
        return BernsteinPolynomial(
            check_ptr(gelfgren_bernstein_subtract(ptr_.get(), other.ptr_.get())));
    }

    /**
     * @brief Multiply two polynomials
     */
    BernsteinPolynomial operator*(const BernsteinPolynomial& other) const {
        return BernsteinPolynomial(
            check_ptr(gelfgren_bernstein_multiply(ptr_.get(), other.ptr_.get())));
    }

    /**
     * @brief Get raw C pointer (for advanced usage)
     */
    GelfgrenBernstein* get() const { return ptr_.get(); }

private:
    // Private constructor for wrapping existing C pointers
    explicit BernsteinPolynomial(GelfgrenBernstein* ptr)
        : ptr_(check_ptr(ptr), &gelfgren_bernstein_free) {}

    std::unique_ptr<GelfgrenBernstein, decltype(&gelfgren_bernstein_free)> ptr_;

    friend class RationalFunction;
};

/**
 * @brief RAII wrapper for rational function P(x)/Q(x)
 *
 * Represents a rational function as ratio of two Bernstein polynomials.
 *
 * @example
 * ```cpp
 * BernsteinPolynomial num({1.0, 2.0}, 0.0, 1.0);
 * BernsteinPolynomial den({1.0, 1.0}, 0.0, 1.0);
 * RationalFunction r(num, den);
 * double value = r(0.5);  // May throw if denominator is zero
 * ```
 */
class RationalFunction {
public:
    /**
     * @brief Construct from numerator and denominator polynomials
     */
    RationalFunction(const BernsteinPolynomial& numerator,
                     const BernsteinPolynomial& denominator)
        : ptr_(check_ptr(gelfgren_rational_create(
            numerator.get(), denominator.get())),
            &gelfgren_rational_free) {}

    /**
     * @brief Evaluate rational function at point x
     * @throws GelfgrenException if denominator is zero at x
     */
    double operator()(double x) const {
        double result;
        check_error(gelfgren_rational_evaluate(ptr_.get(), x, &result));
        return result;
    }

    /**
     * @brief Evaluate rational function (named method)
     */
    double evaluate(double x) const {
        return (*this)(x);
    }

    /**
     * @brief Compute derivative using quotient rule
     */
    RationalFunction derivative() const {
        return RationalFunction(
            check_ptr(gelfgren_rational_derivative(ptr_.get())));
    }

    /**
     * @brief Get raw C pointer (for advanced usage)
     */
    GelfgrenRational* get() const { return ptr_.get(); }

private:
    explicit RationalFunction(GelfgrenRational* ptr)
        : ptr_(check_ptr(ptr), &gelfgren_rational_free) {}

    std::unique_ptr<GelfgrenRational, decltype(&gelfgren_rational_free)> ptr_;
};

/**
 * @brief RAII wrapper for Padé approximant [n/m]
 *
 * Padé approximants provide rational approximations to functions from
 * power series or derivative data.
 *
 * @example
 * ```cpp
 * // Approximate exp(x) with [2/2] Padé
 * std::vector<double> exp_series = {1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0};
 * PadeApproximant pade(exp_series, 2, 2, 0.0, -1.0, 1.0);
 * double approx = pade(0.5);  // ≈ exp(0.5)
 * ```
 */
class PadeApproximant {
public:
    /**
     * @brief Construct Padé approximant from power series
     * @param coeffs Power series coefficients [c₀, c₁, ..., c_{n+m}]
     * @param n Numerator degree
     * @param m Denominator degree
     * @param center Expansion center
     * @param a Left endpoint of interval
     * @param b Right endpoint of interval
     */
    PadeApproximant(const std::vector<double>& coeffs,
                    size_t n, size_t m, double center, double a, double b)
        : ptr_(check_ptr(gelfgren_pade_from_series(
            coeffs.data(), n, m, center, a, b)),
            &gelfgren_pade_free) {}

    /**
     * @brief Evaluate Padé approximant at point x
     */
    double operator()(double x) const {
        double result;
        check_error(gelfgren_pade_evaluate(ptr_.get(), x, &result));
        return result;
    }

    /**
     * @brief Evaluate Padé approximant (named method)
     */
    double evaluate(double x) const {
        return (*this)(x);
    }

private:
    std::unique_ptr<GelfgrenPade, decltype(&gelfgren_pade_free)> ptr_;
};

/**
 * @brief RAII wrapper for mesh partition
 *
 * Represents a partition of an interval [a, b] for piecewise approximation.
 */
class Mesh {
public:
    /**
     * @brief Create uniform mesh with equal spacing
     * @param a Left endpoint
     * @param b Right endpoint
     * @param n Number of subintervals
     */
    static Mesh uniform(double a, double b, size_t n) {
        return Mesh(check_ptr(gelfgren_mesh_uniform(a, b, n)));
    }

    /**
     * @brief Create Chebyshev mesh with boundary clustering
     * @param a Left endpoint
     * @param b Right endpoint
     * @param n Number of subintervals
     */
    static Mesh chebyshev(double a, double b, size_t n) {
        return Mesh(check_ptr(gelfgren_mesh_chebyshev(a, b, n)));
    }

    /**
     * @brief Get raw C pointer (for advanced usage)
     */
    GelfgrenMesh* get() const { return ptr_.get(); }

private:
    explicit Mesh(GelfgrenMesh* ptr)
        : ptr_(check_ptr(ptr), &gelfgren_mesh_free) {}

    std::unique_ptr<GelfgrenMesh, decltype(&gelfgren_mesh_free)> ptr_;
};

} // namespace gelfgren
