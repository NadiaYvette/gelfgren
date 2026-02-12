(* Basic usage example for Gelfgren OCaml bindings *)

open Gelfgren

let () =
  Printf.printf "Gelfgren OCaml Bindings Example\n";
  Printf.printf "%s\n\n" (String.make 60 '=');

  (* Example 1: Bernstein Polynomial *)
  Printf.printf "1. Bernstein Polynomial\n";
  Printf.printf "%s\n" (String.make 60 '-');

  (* Using with_poly for automatic cleanup *)
  BernsteinPolynomial.with_poly [|1.0; 2.0; 6.0|] 0.0 1.0 (fun poly ->
      (* Get degree *)
      let deg = BernsteinPolynomial.degree poly in
      Printf.printf "Created polynomial with degree: %d\n" deg;

      (* Evaluate at a point *)
      let y = BernsteinPolynomial.evaluate poly 0.5 in
      Printf.printf "f(0.5) = %.6f\n" y;

      (* Vectorized evaluation *)
      Printf.printf "\nVectorized evaluation:\n";
      let xs = [|0.0; 0.25; 0.5; 0.75; 1.0|] in
      let ys = BernsteinPolynomial.evaluate_array poly xs in
      Array.iter2
        (fun x y -> Printf.printf "  f(%.2f) = %.6f\n" x y)
        xs ys;

      (* Derivative *)
      Printf.printf "\n";
      let dpoly = BernsteinPolynomial.derivative poly in
      let dy = BernsteinPolynomial.evaluate dpoly 0.5 in
      Printf.printf "Derivative at x=0.5: %.6f\n" dy;
      BernsteinPolynomial.free dpoly;

      (* Integral *)
      let ipoly = BernsteinPolynomial.integral poly in
      let y1 = BernsteinPolynomial.evaluate ipoly 1.0 in
      let y0 = BernsteinPolynomial.evaluate ipoly 0.0 in
      let area = y1 -. y0 in
      Printf.printf "Integral from 0 to 1: %.6f\n" area;
      BernsteinPolynomial.free ipoly
    );

  Printf.printf "\n";

  (* Example 2: Polynomial Arithmetic *)
  Printf.printf "2. Polynomial Arithmetic\n";
  Printf.printf "%s\n" (String.make 60 '-');

  let p1 = BernsteinPolynomial.create [|1.0; 2.0; 3.0|] 0.0 1.0 in
  let p2 = BernsteinPolynomial.create [|2.0; 1.0; 1.0|] 0.0 1.0 in

  let y1 = BernsteinPolynomial.evaluate p1 0.5 in
  let y2 = BernsteinPolynomial.evaluate p2 0.5 in
  Printf.printf "p1(0.5) = %.6f\n" y1;
  Printf.printf "p2(0.5) = %.6f\n" y2;

  (* Addition using infix operator *)
  let open BernsteinPolynomial in
  let sum = p1 + p2 in
  let y_sum = evaluate sum 0.5 in
  Printf.printf "p1 + p2 at x=0.5: %.6f\n" y_sum;
  free sum;

  (* Scaling using infix operator *)
  let scaled = 2.0 *. p1 in
  let y_scaled = evaluate scaled 0.5 in
  Printf.printf "2 * p1 at x=0.5: %.6f\n" y_scaled;
  free scaled;

  free p1;
  free p2;

  Printf.printf "\n";

  (* Example 3: Rational Function *)
  Printf.printf "3. Rational Function\n";
  Printf.printf "%s\n" (String.make 60 '-');

  (* Create rational function: x / (1 + x) on [0, 1] *)
  let num = BernsteinPolynomial.create [|0.0; 1.0|] 0.0 1.0 in
  let den = BernsteinPolynomial.create [|1.0; 2.0|] 0.0 1.0 in

  RationalFunction.with_rational num den (fun rat ->
      let x = 0.5 in
      let y = RationalFunction.evaluate rat x in
      let expected = x /. (1.0 +. x) in
      Printf.printf "R(%.1f) = %.6f\n" x y;
      Printf.printf "Expected: %.6f\n" expected;

      (* Vectorized evaluation *)
      Printf.printf "\nVectorized evaluation:\n";
      let xs = [|0.0; 0.25; 0.5; 0.75; 1.0|] in
      let ys = RationalFunction.evaluate_array rat xs in
      Array.iter2 (fun x y ->
          let expected = x /. (1.0 +. x) in
          Printf.printf "  R(%.2f) = %.6f (expected: %.6f)\n" x y expected
        ) xs ys
    );

  BernsteinPolynomial.free num;
  BernsteinPolynomial.free den;

  Printf.printf "\n";

  (* Example 4: Padé Approximant *)
  Printf.printf "4. Padé Approximant\n";
  Printf.printf "%s\n" (String.make 60 '-');

  (* Approximate exp(x) near x=0 with [2/2] Padé approximant *)
  let coeffs = [|1.0; 1.0; 0.5; 1.0/.6.0; 1.0/.24.0|] in

  PadeApproximant.with_pade coeffs 2 2 (-1.0) 1.0 (fun pade ->
      Printf.printf "Approximation of exp(x):\n";
      let xs = [|-0.5; -0.25; 0.0; 0.25; 0.5|] in
      let ys = PadeApproximant.evaluate_array pade xs in
      Array.iter2 (fun x y ->
          let exact = exp x in
          let error = abs_float (y -. exact) in
          Printf.printf "  x=%.2f: Padé=%.6f, exp(x)=%.6f, error=%.8f\n"
            x y exact error
        ) xs ys
    );

  Printf.printf "\n";
  Printf.printf "All examples completed successfully!\n"
