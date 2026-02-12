(* OCaml Ctypes bindings for Gelfgren C FFI *)

open Ctypes
open Foreign

(* Error codes *)
type error_code =
  | Success
  | NullPointer
  | InvalidInterval
  | EmptyCoefficients
  | SingularMatrix
  | DivisionByZero
  | Pole
  | Unknown

let error_code_of_int = function
  | 0 -> Success
  | -1 -> NullPointer
  | -2 -> InvalidInterval
  | -3 -> EmptyCoefficients
  | -4 -> SingularMatrix
  | -5 -> DivisionByZero
  | -6 -> Pole
  | _ -> Unknown

(* Opaque pointer types *)
type bernstein_polynomial
let bernstein_polynomial : bernstein_polynomial structure typ = structure "GelfgrenBernstein"

type rational_function
let rational_function : rational_function structure typ = structure "GelfgrenRational"

type pade_approximant
let pade_approximant : pade_approximant structure typ = structure "GelfgrenPade"

(* Error handling *)
let last_error_message =
  foreign "gelfgren_last_error_message"
    (void @-> returning string)

(* Bernstein polynomial functions *)
let bernstein_create =
  foreign "gelfgren_bernstein_create"
    (ptr double @-> size_t @-> double @-> double @->
     returning (ptr bernstein_polynomial))

let bernstein_free =
  foreign "gelfgren_bernstein_free"
    (ptr bernstein_polynomial @-> returning void)

let bernstein_evaluate =
  foreign "gelfgren_bernstein_evaluate"
    (ptr bernstein_polynomial @-> double @-> ptr double @->
     returning int)

let bernstein_derivative =
  foreign "gelfgren_bernstein_derivative"
    (ptr bernstein_polynomial @-> returning (ptr bernstein_polynomial))

let bernstein_integral =
  foreign "gelfgren_bernstein_integral"
    (ptr bernstein_polynomial @-> returning (ptr bernstein_polynomial))

let bernstein_degree =
  foreign "gelfgren_bernstein_degree"
    (ptr bernstein_polynomial @-> ptr int @-> returning int)

let bernstein_add =
  foreign "gelfgren_bernstein_add"
    (ptr bernstein_polynomial @-> ptr bernstein_polynomial @->
     returning (ptr bernstein_polynomial))

let bernstein_subtract =
  foreign "gelfgren_bernstein_subtract"
    (ptr bernstein_polynomial @-> ptr bernstein_polynomial @->
     returning (ptr bernstein_polynomial))

let bernstein_multiply =
  foreign "gelfgren_bernstein_multiply"
    (ptr bernstein_polynomial @-> ptr bernstein_polynomial @->
     returning (ptr bernstein_polynomial))

let bernstein_scale =
  foreign "gelfgren_bernstein_scale"
    (double @-> ptr bernstein_polynomial @->
     returning (ptr bernstein_polynomial))

(* Rational function functions *)
let rational_create =
  foreign "gelfgren_rational_create"
    (ptr bernstein_polynomial @-> ptr bernstein_polynomial @->
     returning (ptr rational_function))

let rational_free =
  foreign "gelfgren_rational_free"
    (ptr rational_function @-> returning void)

let rational_evaluate =
  foreign "gelfgren_rational_evaluate"
    (ptr rational_function @-> double @-> ptr double @->
     returning int)

let rational_derivative =
  foreign "gelfgren_rational_derivative"
    (ptr rational_function @-> returning (ptr rational_function))

(* Padé approximant functions *)
let pade_from_series =
  foreign "gelfgren_pade_from_series"
    (ptr double @-> size_t @-> int @-> int @-> double @-> double @->
     returning (ptr pade_approximant))

let pade_free =
  foreign "gelfgren_pade_free"
    (ptr pade_approximant @-> returning void)

let pade_evaluate =
  foreign "gelfgren_pade_evaluate"
    (ptr pade_approximant @-> double @-> ptr double @->
     returning int)
