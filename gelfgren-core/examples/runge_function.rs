/// Demonstration of rational approximation for Runge's function
///
/// Runge's function 1/(1+25x²) is a classic example where polynomial
/// interpolation fails catastrophically (Runge's phenomenon) but rational
/// approximation excels.

use gelfgren_core::pade::TwoPointPade;

fn runge(x: f64) -> f64 {
    1.0 / (1.0 + 25.0 * x * x)
}

fn runge_derivative(x: f64) -> f64 {
    -50.0 * x / (1.0 + 25.0 * x * x).powi(2)
}

fn main() {
    println!("Runge's Function: 1/(1 + 25x²)");
    println!("{}", "=".repeat(60));
    println!();

    // Approximate on [0, 1]
    let x0 = 0.0;
    let x1 = 1.0;

    // For [2/1]: n+m+1 = 4, so p=2 (need 2 derivatives at each endpoint)
    let left_derivs = vec![
        runge(x0),
        runge_derivative(x0),
    ];

    let right_derivs = vec![
        runge(x1),
        runge_derivative(x1),
    ];

    println!("Endpoint values:");
    println!("  f({}) = {:.8}, f'({}) = {:.8}", x0, left_derivs[0], x0, left_derivs[1]);
    println!("  f({}) = {:.8}, f'({}) = {:.8}", x1, right_derivs[0], x1, right_derivs[1]);
    println!();

    // Construct [2/1] rational
    println!("Constructing [2/1] rational approximant...");
    let rat_pade = TwoPointPade::from_endpoint_derivatives(
        &left_derivs,
        &right_derivs,
        2,  // numerator degree
        1,  // denominator degree
        x0,
        x1,
    ).unwrap();

    println!();
    println!("{:>6} {:>12} {:>12} {:>12}",
             "x", "Exact", "[2/1] Rat", "Error");
    println!("{}", "-".repeat(48));

    // Evaluate at test points
    let test_points: Vec<f64> = (0..=20).map(|i| i as f64 * 0.05).collect();

    let mut max_err: f64 = 0.0;
    let mut sum_err2 = 0.0;

    for &x in &test_points {
        let exact = runge(x);
        let rat_approx = rat_pade.evaluate(x).unwrap();
        let err = (rat_approx - exact).abs();

        max_err = max_err.max(err);
        sum_err2 += err * err;

        println!("{:6.3} {:12.8} {:12.8} {:12.2e}",
                 x, exact, rat_approx, err);
    }

    println!("{}", "-".repeat(48));

    let l2_err = (sum_err2 / test_points.len() as f64).sqrt();

    println!();
    println!("Error Summary:");
    println!("  L∞ error: {:.6e}", max_err);
    println!("  L2 error: {:.6e}", l2_err);
    println!();

    println!("Denominator coefficients (Bernstein basis):");
    let denom = rat_pade.rational().denominator();
    let denom_coeffs = denom.scaled_coefficients();
    println!("  {:?}", denom_coeffs);

    // Check if denominator is non-trivial
    let is_constant = denom_coeffs.windows(2).all(|w| (w[0] - w[1]).abs() < 1e-10);
    if is_constant {
        println!("  WARNING: Denominator is constant!");
    } else {
        let variation = denom_coeffs.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap()
            - denom_coeffs.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        println!("  ✓ Denominator varies by {:.6e}", variation);
    }

    println!();
    println!("Numerator coefficients (Bernstein basis):");
    let numer = rat_pade.rational().numerator();
    let numer_coeffs = numer.scaled_coefficients();
    println!("  {:?}", numer_coeffs);
}
