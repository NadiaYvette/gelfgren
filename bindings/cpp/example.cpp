/**
 * @file example.cpp
 * @brief Example usage of Gelfgren C++ API
 *
 * Demonstrates object-oriented interface with RAII and operator overloading.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "gelfgren.hpp"

using namespace gelfgren;

void demonstrate_polynomials() {
    std::cout << "Bernstein Polynomial Operations\n";
    std::cout << "================================\n\n";

    // Create polynomial P(x) = 1 + 2x + 3x²
    BernsteinPolynomial p({1.0, 2.0, 3.0}, 0.0, 1.0);

    std::cout << "Created polynomial P(x) of degree " << p.degree() << "\n";
    auto [a, b] = p.interval();
    std::cout << "Interval: [" << a << ", " << b << "]\n\n";

    // Evaluate at several points
    std::cout << "Evaluations:\n";
    for (double x : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        std::cout << "  P(" << std::fixed << std::setprecision(2) << x
                  << ") = " << std::setprecision(6) << p(x) << "\n";
    }

    // Compute derivative using member function
    std::cout << "\nDerivative P'(x):\n";
    auto p_prime = p.derivative();
    for (double x : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        std::cout << "  P'(" << std::fixed << std::setprecision(2) << x
                  << ") = " << std::setprecision(6) << p_prime(x) << "\n";
    }

    // Demonstrate polynomial arithmetic
    std::cout << "\nPolynomial Arithmetic:\n";
    BernsteinPolynomial q({1.0, 0.0}, 0.0, 1.0);  // q(x) = 1
    auto sum = p + q;
    auto diff = p - q;
    auto prod = p * q;

    double x = 0.5;
    std::cout << "  At x = " << x << ":\n";
    std::cout << "    P(x) + Q(x) = " << sum(x) << "\n";
    std::cout << "    P(x) - Q(x) = " << diff(x) << "\n";
    std::cout << "    P(x) * Q(x) = " << prod(x) << "\n";
}

void demonstrate_rational_functions() {
    std::cout << "\n\nRational Functions\n";
    std::cout << "==================\n\n";

    // Create R(x) = (1 + x) / (1 + 2x)
    BernsteinPolynomial num({1.0, 1.0}, 0.0, 1.0);   // 1 + x
    BernsteinPolynomial den({1.0, 2.0}, 0.0, 1.0);   // 1 + 2x

    RationalFunction r(num, den);

    std::cout << "Created rational function R(x) = (1 + x) / (1 + 2x)\n\n";

    std::cout << "Evaluations:\n";
    for (double x : {0.0, 0.2, 0.4, 0.6, 0.8, 1.0}) {
        std::cout << "  R(" << std::fixed << std::setprecision(1) << x
                  << ") = " << std::setprecision(6) << r(x) << "\n";
    }

    // Derivative using quotient rule
    std::cout << "\nDerivative R'(x):\n";
    auto r_prime = r.derivative();
    for (double x : {0.0, 0.2, 0.4, 0.6, 0.8, 1.0}) {
        std::cout << "  R'(" << std::fixed << std::setprecision(1) << x
                  << ") = " << std::setprecision(6) << r_prime(x) << "\n";
    }
}

void demonstrate_pade_approximants() {
    std::cout << "\n\nPadé Approximants\n";
    std::cout << "=================\n\n";

    // Approximate exp(x) around x = 0 using Taylor series
    // exp(x) = 1 + x + x²/2 + x³/6 + x⁴/24 + ...
    std::vector<double> exp_coeffs = {
        1.0,           // x⁰
        1.0,           // x¹
        0.5,           // x²/2!
        1.0/6.0,       // x³/3!
        1.0/24.0       // x⁴/4!
    };

    // Create [2/2] Padé approximant
    PadeApproximant pade(exp_coeffs, 2, 2, 0.0, -1.0, 1.0);

    std::cout << "Padé [2/2] approximant for exp(x) on [-1, 1]\n\n";
    std::cout << std::setw(8) << "x"
              << std::setw(15) << "exp(x)"
              << std::setw(15) << "Padé(x)"
              << std::setw(15) << "Error\n";
    std::cout << std::string(52, '-') << "\n";

    for (double x = -1.0; x <= 1.0; x += 0.25) {
        double exact = std::exp(x);
        double approx = pade(x);
        double error = std::abs(exact - approx);

        std::cout << std::fixed << std::setprecision(2) << std::setw(8) << x
                  << std::setprecision(6) << std::setw(15) << exact
                  << std::setw(15) << approx
                  << std::scientific << std::setw(15) << error << "\n";
    }
}

void demonstrate_meshes() {
    std::cout << "\n\nMesh Generation\n";
    std::cout << "===============\n\n";

    // Create uniform mesh
    auto uniform_mesh = Mesh::uniform(0.0, 1.0, 4);
    std::cout << "Created uniform mesh on [0, 1] with 4 subintervals\n";

    // Create Chebyshev mesh
    auto cheb_mesh = Mesh::chebyshev(0.0, 1.0, 4);
    std::cout << "Created Chebyshev mesh on [0, 1] with 4 subintervals\n";
    std::cout << "  (Chebyshev meshes cluster points near boundaries)\n";
}

void demonstrate_exception_handling() {
    std::cout << "\n\nException Handling\n";
    std::cout << "==================\n\n";

    try {
        // Create a rational function
        BernsteinPolynomial num({1.0, 1.0}, 0.0, 1.0);
        BernsteinPolynomial den({1.0, -1.0}, 0.0, 1.0);  // Zero at x = 0.5
        RationalFunction r(num, den);

        std::cout << "Evaluating rational function with pole...\n";

        // This should work
        std::cout << "  R(0.0) = " << r(0.0) << " (OK)\n";

        // This should throw (pole at x ≈ 0.5)
        std::cout << "  R(0.5) = ";
        double result = r(0.5);
        std::cout << result << "\n";

    } catch (const GelfgrenException& e) {
        std::cout << "  Caught exception: " << e.what() << "\n";
        std::cout << "  (This is expected behavior for poles)\n";
    }
}

int main() {
    std::cout << "Gelfgren C++ API Example\n";
    std::cout << "========================\n\n";

    try {
        demonstrate_polynomials();
        demonstrate_rational_functions();
        demonstrate_pade_approximants();
        demonstrate_meshes();
        demonstrate_exception_handling();

        std::cout << "\n\nAll examples completed successfully!\n";
        return 0;

    } catch (const GelfgrenException& e) {
        std::cerr << "\nError: " << e.what() << "\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "\nUnexpected error: " << e.what() << "\n";
        return 1;
    }
}
