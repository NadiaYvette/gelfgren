(* Unit tests for Gelfgren OCaml bindings *)

open OUnit2
open Gelfgren

let epsilon = 1e-10

let assert_float_equal = assert_equal ~cmp:(fun a b -> abs_float (a -. b) < epsilon) ~printer:string_of_float

(* Test Bernstein Polynomial *)
let test_bernstein_create _ctx =
  let poly = BernsteinPolynomial.create [|1.0; 2.0; 6.0|] 0.0 1.0 in
  let deg = BernsteinPolynomial.degree poly in
  assert_equal 2 deg;
  BernsteinPolynomial.free poly

let test_bernstein_evaluate _ctx =
  BernsteinPolynomial.with_poly [|1.0; 2.0; 6.0|] 0.0 1.0 (fun poly ->
      let y = BernsteinPolynomial.evaluate poly 0.5 in
      assert_float_equal 3.0 y
    )

let test_bernstein_derivative _ctx =
  BernsteinPolynomial.with_poly [|1.0; 2.0; 6.0|] 0.0 1.0 (fun poly ->
      let dpoly = BernsteinPolynomial.derivative poly in
      let y = BernsteinPolynomial.evaluate dpoly 0.5 in
      (* Derivative of 2x^2 + 3x + 1 is 4x + 3, at x=0.5 is 5.0 *)
      assert_float_equal 5.0 y;
      BernsteinPolynomial.free dpoly
    )

let test_bernstein_integral _ctx =
  BernsteinPolynomial.with_poly [|1.0; 2.0; 6.0|] 0.0 1.0 (fun poly ->
      let ipoly = BernsteinPolynomial.integral poly in
      let y1 = BernsteinPolynomial.evaluate ipoly 1.0 in
      let y0 = BernsteinPolynomial.evaluate ipoly 0.0 in
      let area = y1 -. y0 in
      (* Integral of 2x^2 + 3x + 1 from 0 to 1 is [2/3 x^3 + 3/2 x^2 + x] = 2/3 + 3/2 + 1 *)
      let expected = (2.0 /. 3.0) +. (3.0 /. 2.0) +. 1.0 in
      assert_float_equal expected area;
      BernsteinPolynomial.free ipoly
    )

let test_bernstein_add _ctx =
  let p1 = BernsteinPolynomial.create [|1.0; 2.0; 3.0|] 0.0 1.0 in
  let p2 = BernsteinPolynomial.create [|2.0; 1.0; 1.0|] 0.0 1.0 in
  let open BernsteinPolynomial in
  let sum = p1 + p2 in
  let y = evaluate sum 0.5 in
  let y1 = evaluate p1 0.5 in
  let y2 = evaluate p2 0.5 in
  assert_float_equal (y1 +. y2) y;
  free p1;
  free p2;
  free sum

let test_bernstein_scale _ctx =
  let p = BernsteinPolynomial.create [|1.0; 2.0; 3.0|] 0.0 1.0 in
  let open BernsteinPolynomial in
  let scaled = 2.0 *. p in
  let y_orig = evaluate p 0.5 in
  let y_scaled = evaluate scaled 0.5 in
  assert_float_equal (2.0 *. y_orig) y_scaled;
  free p;
  free scaled

let test_bernstein_multiply _ctx =
  let p1 = BernsteinPolynomial.create [|1.0; 2.0|] 0.0 1.0 in
  let p2 = BernsteinPolynomial.create [|1.0; 1.0|] 0.0 1.0 in
  let open BernsteinPolynomial in
  let prod = p1 * p2 in
  let y = evaluate prod 0.5 in
  let y1 = evaluate p1 0.5 in
  let y2 = evaluate p2 0.5 in
  assert_float_equal (y1 *. y2) y;
  free p1;
  free p2;
  free prod

(* Test Rational Functions *)
let test_rational_create _ctx =
  let num = BernsteinPolynomial.create [|0.0; 1.0|] 0.0 1.0 in
  let den = BernsteinPolynomial.create [|1.0; 2.0|] 0.0 1.0 in
  let rat = RationalFunction.create num den in
  RationalFunction.free rat;
  BernsteinPolynomial.free num;
  BernsteinPolynomial.free den

let test_rational_evaluate _ctx =
  let num = BernsteinPolynomial.create [|0.0; 1.0|] 0.0 1.0 in
  let den = BernsteinPolynomial.create [|1.0; 2.0|] 0.0 1.0 in
  RationalFunction.with_rational num den (fun rat ->
      let x = 0.5 in
      let y = RationalFunction.evaluate rat x in
      let expected = x /. (1.0 +. x) in
      assert_float_equal expected y
    );
  BernsteinPolynomial.free num;
  BernsteinPolynomial.free den

let test_rational_array _ctx =
  let num = BernsteinPolynomial.create [|0.0; 1.0|] 0.0 1.0 in
  let den = BernsteinPolynomial.create [|1.0; 2.0|] 0.0 1.0 in
  RationalFunction.with_rational num den (fun rat ->
      let xs = [|0.0; 0.5; 1.0|] in
      let ys = RationalFunction.evaluate_array rat xs in
      Array.iteri (fun i x ->
          let expected = x /. (1.0 +. x) in
          assert_float_equal expected ys.(i)
        ) xs
    );
  BernsteinPolynomial.free num;
  BernsteinPolynomial.free den

(* Test Padé Approximants *)
let test_pade_from_series _ctx =
  let coeffs = [|1.0; 1.0; 0.5; 1.0/.6.0; 1.0/.24.0|] in
  let pade = PadeApproximant.from_series coeffs 2 2 (-1.0) 1.0 in
  PadeApproximant.free pade

let test_pade_evaluate _ctx =
  let coeffs = [|1.0; 1.0; 0.5; 1.0/.6.0; 1.0/.24.0|] in
  PadeApproximant.with_pade coeffs 2 2 (-1.0) 1.0 (fun pade ->
      let x = 0.5 in
      let y = PadeApproximant.evaluate pade x in
      let exact = exp x in
      let error = abs_float (y -. exact) in
      (* Padé [2/2] should be quite accurate for exp near 0 *)
      assert_bool "Padé error too large" (error < 0.001)
    )

let test_pade_array _ctx =
  let coeffs = [|1.0; 1.0; 0.5; 1.0/.6.0; 1.0/.24.0|] in
  PadeApproximant.with_pade coeffs 2 2 (-1.0) 1.0 (fun pade ->
      let xs = [|-0.5; 0.0; 0.5|] in
      let ys = PadeApproximant.evaluate_array pade xs in
      Array.iteri (fun i x ->
          let exact = exp x in
          let error = abs_float (ys.(i) -. exact) in
          assert_bool "Padé error too large" (error < 0.001)
        ) xs
    )

(* Test error handling *)
let test_empty_coefficients _ctx =
  assert_raises (Gelfgren_error "Empty coefficient array")
    (fun () -> BernsteinPolynomial.create [||] 0.0 1.0)

let test_freed_polynomial _ctx =
  let poly = BernsteinPolynomial.create [|1.0; 2.0|] 0.0 1.0 in
  BernsteinPolynomial.free poly;
  assert_raises (Gelfgren_error "Polynomial already freed")
    (fun () -> BernsteinPolynomial.evaluate poly 0.5 |> ignore)

(* Test suite *)
let suite =
  "Gelfgren Tests" >::: [
    "test_bernstein_create" >:: test_bernstein_create;
    "test_bernstein_evaluate" >:: test_bernstein_evaluate;
    "test_bernstein_derivative" >:: test_bernstein_derivative;
    "test_bernstein_integral" >:: test_bernstein_integral;
    "test_bernstein_add" >:: test_bernstein_add;
    "test_bernstein_scale" >:: test_bernstein_scale;
    "test_bernstein_multiply" >:: test_bernstein_multiply;
    "test_rational_create" >:: test_rational_create;
    "test_rational_evaluate" >:: test_rational_evaluate;
    "test_rational_array" >:: test_rational_array;
    "test_pade_from_series" >:: test_pade_from_series;
    "test_pade_evaluate" >:: test_pade_evaluate;
    "test_pade_array" >:: test_pade_array;
    "test_empty_coefficients" >:: test_empty_coefficients;
    "test_freed_polynomial" >:: test_freed_polynomial;
  ]

let () = run_test_tt_main suite
