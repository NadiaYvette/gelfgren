/**
 * Basic usage example of Gelfgren C API
 *
 * Demonstrates creating and evaluating Bernstein polynomials.
 */

#include <stdio.h>
#include <stdlib.h>
#include "../../include/gelfgren.h"

int main() {
    printf("Gelfgren C API Example\n");
    printf("======================\n\n");

    // Create a Bernstein polynomial: P(x) = 1 + 2x + 3x²
    // Coefficients in Bernstein form: [1, 2, 3]
    double coeffs[] = {1.0, 2.0, 3.0};
    GelfgrenBernstein* poly = gelfgren_bernstein_create(coeffs, 2, 0.0, 1.0);

    if (poly == NULL) {
        const char* error = gelfgren_last_error_message();
        fprintf(stderr, "Error creating polynomial: %s\n", error);
        return 1;
    }

    printf("Created Bernstein polynomial of degree 2 on [0, 1]\n");
    printf("Coefficients: [%.1f, %.1f, %.1f]\n\n", coeffs[0], coeffs[1], coeffs[2]);

    // Get polynomial degree
    uintptr_t degree;
    GelfgrenErrorCode code = gelfgren_bernstein_degree(poly, &degree);
    if (code != SUCCESS) {
        fprintf(stderr, "Error getting degree: %s\n", gelfgren_last_error_message());
        gelfgren_bernstein_free(poly);
        return 1;
    }
    printf("Degree: %zu\n", degree);

    // Get polynomial interval
    double a, b;
    code = gelfgren_bernstein_interval(poly, &a, &b);
    if (code != SUCCESS) {
        fprintf(stderr, "Error getting interval: %s\n", gelfgren_last_error_message());
        gelfgren_bernstein_free(poly);
        return 1;
    }
    printf("Interval: [%.1f, %.1f]\n\n", a, b);

    // Evaluate at several points
    printf("Evaluations:\n");
    double points[] = {0.0, 0.25, 0.5, 0.75, 1.0};
    int num_points = sizeof(points) / sizeof(points[0]);

    for (int i = 0; i < num_points; i++) {
        double x = points[i];
        double result;
        code = gelfgren_bernstein_evaluate(poly, x, &result);

        if (code != SUCCESS) {
            fprintf(stderr, "Error evaluating at x=%.2f: %s\n",
                    x, gelfgren_last_error_message());
            gelfgren_bernstein_free(poly);
            return 1;
        }

        printf("  P(%.2f) = %.6f\n", x, result);
    }

    // Compute derivative
    printf("\nComputing derivative...\n");
    GelfgrenBernstein* deriv = gelfgren_bernstein_derivative(poly);
    if (deriv == NULL) {
        fprintf(stderr, "Error computing derivative: %s\n", gelfgren_last_error_message());
        gelfgren_bernstein_free(poly);
        return 1;
    }

    printf("Derivative evaluations:\n");
    for (int i = 0; i < num_points; i++) {
        double x = points[i];
        double result;
        code = gelfgren_bernstein_evaluate(deriv, x, &result);

        if (code != SUCCESS) {
            fprintf(stderr, "Error evaluating derivative at x=%.2f: %s\n",
                    x, gelfgren_last_error_message());
            gelfgren_bernstein_free(deriv);
            gelfgren_bernstein_free(poly);
            return 1;
        }

        printf("  P'(%.2f) = %.6f\n", x, result);
    }

    // Clean up
    gelfgren_bernstein_free(deriv);
    gelfgren_bernstein_free(poly);
    gelfgren_clear_last_error();

    printf("\nExample completed successfully!\n");
    return 0;
}
