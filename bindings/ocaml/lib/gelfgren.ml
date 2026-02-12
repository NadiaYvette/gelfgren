(* High-level OCaml API for Gelfgren *)

open Ctypes
open Gelfgren_bindings

(* Exception for Gelfgren errors *)
exception Gelfgren_error of string

(* Helper to check error codes and raise exceptions *)
let check_error code =
  match error_code_of_int code with
  | Success -> ()
  | _ ->
      let msg = last_error_message () in
      raise (Gelfgren_error msg)

(* Helper to check null pointers *)
let check_null ptr name =
  if is_null ptr then
    let msg = last_error_message () in
    raise (Gelfgren_error (name ^ ": " ^ msg))
  else ptr

(* Bernstein Polynomial module *)
module BernsteinPolynomial = struct
  type t = {
    ptr : bernstein_polynomial Ctypes.structure Ctypes.ptr;
    mutable freed : bool;
  }

  let create coeffs a b =
    let n = Array.length coeffs in
    if n = 0 then
      raise (Gelfgren_error "Empty coefficient array");
    let degree = n - 1 in
    let coeffs_ptr = CArray.of_list double (Array.to_list coeffs) in
    let ptr = bernstein_create
        (CArray.start coeffs_ptr)
        (Unsigned.Size_t.of_int degree)
        a b
    in
    let ptr = check_null ptr "create_bernstein" in
    { ptr; freed = false }

  let free t =
    if not t.freed then begin
      bernstein_free t.ptr;
      t.freed <- true
    end

  let with_poly coeffs a b f =
    let poly = create coeffs a b in
    Fun.protect ~finally:(fun () -> free poly) (fun () -> f poly)

  let evaluate t x =
    if t.freed then raise (Gelfgren_error "Polynomial already freed");
    let result = allocate double 0.0 in
    let code = bernstein_evaluate t.ptr x result in
    check_error code;
    !@result

  let evaluate_array t xs =
    Array.map (evaluate t) xs

  let derivative t =
    if t.freed then raise (Gelfgren_error "Polynomial already freed");
    let ptr = bernstein_derivative t.ptr in
    let ptr = check_null ptr "derivative" in
    { ptr; freed = false }

  let integral t =
    if t.freed then raise (Gelfgren_error "Polynomial already freed");
    let ptr = bernstein_integral t.ptr in
    let ptr = check_null ptr "integral" in
    { ptr; freed = false }

  let degree t =
    if t.freed then raise (Gelfgren_error "Polynomial already freed");
    let deg = allocate int 0 in
    let code = bernstein_degree t.ptr deg in
    check_error code;
    !@deg

  let add t1 t2 =
    if t1.freed || t2.freed then
      raise (Gelfgren_error "Polynomial already freed");
    let ptr = bernstein_add t1.ptr t2.ptr in
    let ptr = check_null ptr "add" in
    { ptr; freed = false }

  let subtract t1 t2 =
    if t1.freed || t2.freed then
      raise (Gelfgren_error "Polynomial already freed");
    let ptr = bernstein_subtract t1.ptr t2.ptr in
    let ptr = check_null ptr "subtract" in
    { ptr; freed = false }

  let multiply t1 t2 =
    if t1.freed || t2.freed then
      raise (Gelfgren_error "Polynomial already freed");
    let ptr = bernstein_multiply t1.ptr t2.ptr in
    let ptr = check_null ptr "multiply" in
    { ptr; freed = false }

  let scale scalar t =
    if t.freed then raise (Gelfgren_error "Polynomial already freed");
    let ptr = bernstein_scale scalar t.ptr in
    let ptr = check_null ptr "scale" in
    { ptr; freed = false }

  (* Operator overloading *)
  let ( + ) = add
  let ( - ) = subtract
  let ( * ) = multiply
  let ( *. ) = scale
end

(* Rational Function module *)
module RationalFunction = struct
  type t = {
    ptr : rational_function Ctypes.structure Ctypes.ptr;
    mutable freed : bool;
  }

  let create num den =
    if num.BernsteinPolynomial.freed || den.BernsteinPolynomial.freed then
      raise (Gelfgren_error "Polynomial already freed");
    let ptr = rational_create num.BernsteinPolynomial.ptr den.BernsteinPolynomial.ptr in
    let ptr = check_null ptr "create_rational" in
    { ptr; freed = false }

  let free t =
    if not t.freed then begin
      rational_free t.ptr;
      t.freed <- true
    end

  let with_rational num den f =
    let rat = create num den in
    Fun.protect ~finally:(fun () -> free rat) (fun () -> f rat)

  let evaluate t x =
    if t.freed then raise (Gelfgren_error "Rational function already freed");
    let result = allocate double 0.0 in
    let code = rational_evaluate t.ptr x result in
    check_error code;
    !@result

  let evaluate_array t xs =
    Array.map (evaluate t) xs

  let derivative t =
    if t.freed then raise (Gelfgren_error "Rational function already freed");
    let ptr = rational_derivative t.ptr in
    let ptr = check_null ptr "derivative" in
    { ptr; freed = false }
end

(* Padé Approximant module *)
module PadeApproximant = struct
  type t = {
    ptr : pade_approximant Ctypes.structure Ctypes.ptr;
    mutable freed : bool;
  }

  let from_series coeffs n m a b =
    let len = Array.length coeffs in
    let coeffs_ptr = CArray.of_list double (Array.to_list coeffs) in
    let ptr = pade_from_series
        (CArray.start coeffs_ptr)
        (Unsigned.Size_t.of_int len)
        n m a b
    in
    let ptr = check_null ptr "from_series" in
    { ptr; freed = false }

  let free t =
    if not t.freed then begin
      pade_free t.ptr;
      t.freed <- true
    end

  let with_pade coeffs n m a b f =
    let pade = from_series coeffs n m a b in
    Fun.protect ~finally:(fun () -> free pade) (fun () -> f pade)

  let evaluate t x =
    if t.freed then raise (Gelfgren_error "Padé approximant already freed");
    let result = allocate double 0.0 in
    let code = pade_evaluate t.ptr x result in
    check_error code;
    !@result

  let evaluate_array t xs =
    Array.map (evaluate t) xs
end
