/// Demonstration of rational vs polynomial approximation for e^x
///
/// The exponential function is well-suited for Padé approximation.

use gelfgren_core::pade::TwoPointPade;

fn main() {
    println!("Exponential Function Approximation: e^x on [0, 1]");
    println!("{}", "=".repeat(60));
    println!();

    let x0: f64 = 0.0;
    let x1: f64 = 1.0;

    // For [2/1]: n+m+1 = 4, so p=2 (need 2 derivatives at each endpoint)
    // e^x derivatives are all e^x
    let left_derivs = vec![
        x0.exp(),  // e^0 = 1
        x0.exp(),  // d/dx e^x |_{x=0} = 1
    ];

    let right_derivs = vec![
        x1.exp(),  // e^1 = e
        x1.exp(),  // d/dx e^x |_{x=1} = e
    ];

    println!("Constructing [2/1] rational approximant...");
    let rat_pade = TwoPointPade::from_endpoint_derivatives(
        &left_derivs,
        &right_derivs,
        2,  // numerator degree
        1,  // denominator degree
        x0,
        x1,
    ).unwrap();

    // Also try [3/2] for comparison
    let left_derivs_p3 = vec![
        x0.exp(),
        x0.exp(),
        x0.exp(),
    ];

    let right_derivs_p3 = vec![
        x1.exp(),
        x1.exp(),
        x1.exp(),
    ];

    println!("Constructing [3/2] rational approximant...");
    let rat_pade_32 = TwoPointPade::from_endpoint_derivatives(
        &left_derivs_p3,
        &right_derivs_p3,
        3,  // numerator degree
        2,  // denominator degree
        x0,
        x1,
    ).unwrap();

    println!();
    println!("{:>6} {:>12} {:>12} {:>12} {:>12} {:>12}",
             "x", "Exact", "[2/1] Rat", "[3/2] Rat", "Err [2/1]", "Err [3/2]");
    println!("{}", "-".repeat(78));

    let test_points: Vec<f64> = (0..=20).map(|i| i as f64 * 0.05).collect();

    let mut max_err_21: f64 = 0.0;
    let mut max_err_32: f64 = 0.0;
    let mut sum_err2_21 = 0.0;
    let mut sum_err2_32 = 0.0;

    for &x in &test_points {
        let exact = x.exp();
        let approx_21 = rat_pade.evaluate(x).unwrap();
        let approx_32 = rat_pade_32.evaluate(x).unwrap();

        let err_21 = (approx_21 - exact).abs();
        let err_32 = (approx_32 - exact).abs();

        max_err_21 = max_err_21.max(err_21);
        max_err_32 = max_err_32.max(err_32);
        sum_err2_21 += err_21 * err_21;
        sum_err2_32 += err_32 * err_32;

        println!("{:6.3} {:12.8} {:12.8} {:12.8} {:12.2e} {:12.2e}",
                 x, exact, approx_21, approx_32, err_21, err_32);
    }

    println!("{}", "-".repeat(78));

    let l2_21 = (sum_err2_21 / test_points.len() as f64).sqrt();
    let l2_32 = (sum_err2_32 / test_points.len() as f64).sqrt();

    println!();
    println!("Error Summary:");
    println!("  [2/1] Rational:");
    println!("    L∞ error: {:.6e}", max_err_21);
    println!("    L2 error: {:.6e}", l2_21);
    println!();
    println!("  [3/2] Rational:");
    println!("    L∞ error: {:.6e}", max_err_32);
    println!("    L2 error: {:.6e}", l2_32);
    println!();

    println!("[2/1] Denominator coefficients (Bernstein basis):");
    let denom_21 = rat_pade.rational().denominator();
    let denom_coeffs_21 = denom_21.scaled_coefficients();
    println!("  {:?}", denom_coeffs_21);

    let is_constant = denom_coeffs_21.windows(2).all(|w| (w[0] - w[1]).abs() < 1e-10);
    if is_constant {
        println!("  WARNING: Denominator is constant!");
    } else {
        let variation = denom_coeffs_21.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap()
            - denom_coeffs_21.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        println!("  ✓ Denominator varies by {:.6e}", variation);
    }

    println!();
    println!("[3/2] Denominator coefficients (Bernstein basis):");
    let denom_32 = rat_pade_32.rational().denominator();
    let denom_coeffs_32 = denom_32.scaled_coefficients();
    println!("  {:?}", denom_coeffs_32);

    let is_constant_32 = denom_coeffs_32.windows(2).all(|w| (w[0] - w[1]).abs() < 1e-10);
    if is_constant_32 {
        println!("  WARNING: Denominator is constant!");
    } else {
        let variation = denom_coeffs_32.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap()
            - denom_coeffs_32.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        println!("  ✓ Denominator varies by {:.6e}", variation);
    }
}
