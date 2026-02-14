/// Demonstration of rational vs polynomial Padé approximants
///
/// This example shows that the rational construction produces different
/// (and better) results than polynomial-only approximation.

use gelfgren_core::pade::TwoPointPade;
use std::f64::consts::PI;

fn main() {
    println!("Two-Point Padé Approximation Demonstration");
    println!("{}", "=".repeat(60));
    println!();

    // Approximate f(x) = sin(π x) on [0, 1]
    // We'll compare [1/0] polynomial vs [2/1] rational
    // Note: Must satisfy n+m+1 = 2p (even) for two-point Padé

    let x0 = 0.0;
    let x1 = 1.0;

    // For [1/0]: n+m+1 = 2, so p=1 (need 1 derivative at each endpoint)
    let left_derivs_p1 = vec![
        (PI * x0).sin(),              // f(x0) = sin(0) = 0
    ];

    let right_derivs_p1 = vec![
        (PI * x1).sin(),              // f(x1) = sin(π) = 0
    ];

    // For [2/1]: n+m+1 = 4, so p=2 (need 2 derivatives at each endpoint)
    let left_derivs_p2 = vec![
        (PI * x0).sin(),              // f(x0) = sin(0) = 0
        PI * (PI * x0).cos(),         // f'(x0) = π cos(0) = π
    ];

    let right_derivs_p2 = vec![
        (PI * x1).sin(),              // f(x1) = sin(π) = 0
        PI * (PI * x1).cos(),         // f'(x1) = π cos(π) = -π
    ];

    // Construct [1/0] polynomial (linear)
    println!("Constructing [1/0] polynomial approximant...");
    let poly_pade = TwoPointPade::from_endpoint_derivatives(
        &left_derivs_p1,
        &right_derivs_p1,
        1,  // numerator degree
        0,  // denominator degree (polynomial)
        x0,
        x1,
    ).unwrap();

    // Construct [2/1] rational
    println!("Constructing [2/1] rational approximant...");
    let rat_pade = TwoPointPade::from_endpoint_derivatives(
        &left_derivs_p2,
        &right_derivs_p2,
        2,  // numerator degree
        1,  // denominator degree
        x0,
        x1,
    ).unwrap();

    println!();
    println!("Approximating f(x) = sin(πx) on [0, 1]");
    println!();
    println!("{:>6} {:>12} {:>12} {:>12} {:>12} {:>12}",
             "x", "Exact", "[1/0] Poly", "[2/1] Rat", "Poly Err", "Rat Err");
    println!("{}", "-".repeat(78));

    // Evaluate at test points
    let test_points = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];

    let mut max_poly_err: f64 = 0.0;
    let mut max_rat_err: f64 = 0.0;
    let mut sum_poly_err2 = 0.0;
    let mut sum_rat_err2 = 0.0;

    for &x in &test_points {
        let exact = (PI * x).sin();
        let poly_approx = poly_pade.evaluate(x).unwrap();
        let rat_approx = rat_pade.evaluate(x).unwrap();

        let poly_err = (poly_approx - exact).abs();
        let rat_err = (rat_approx - exact).abs();

        max_poly_err = max_poly_err.max(poly_err);
        max_rat_err = max_rat_err.max(rat_err);
        sum_poly_err2 += poly_err * poly_err;
        sum_rat_err2 += rat_err * rat_err;

        println!("{:6.3} {:12.8} {:12.8} {:12.8} {:12.2e} {:12.2e}",
                 x, exact, poly_approx, rat_approx, poly_err, rat_err);
    }

    println!("{}", "-".repeat(78));

    let l2_poly = (sum_poly_err2 / test_points.len() as f64).sqrt();
    let l2_rat = (sum_rat_err2 / test_points.len() as f64).sqrt();

    println!();
    println!("Error Summary:");
    println!("  Polynomial [1/0]:");
    println!("    L∞ error: {:.6e}", max_poly_err);
    println!("    L2 error: {:.6e}", l2_poly);
    println!();
    println!("  Rational [2/1]:");
    println!("    L∞ error: {:.6e}", max_rat_err);
    println!("    L2 error: {:.6e}", l2_rat);
    println!();

    if max_rat_err < max_poly_err {
        let improvement = ((max_poly_err - max_rat_err) / max_poly_err) * 100.0;
        println!("  Rational improves L∞ error by {:.1}%", improvement);
    } else {
        println!("  Note: Rational and polynomial have similar accuracy for this case");
    }

    println!();
    println!("Denominator coefficients (rational [2/1]):");
    let denom = rat_pade.rational().denominator();
    let denom_coeffs = denom.scaled_coefficients();
    println!("  Bernstein basis: {:?}", denom_coeffs);

    // Check if denominator is non-trivial
    let is_constant = denom_coeffs.windows(2).all(|w| (w[0] - w[1]).abs() < 1e-10);
    if is_constant {
        println!("  WARNING: Denominator is constant! Rational construction may not be working.");
    } else {
        println!("  ✓ Denominator is non-trivial (varying coefficients)");
    }
}
