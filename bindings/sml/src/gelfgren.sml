(* Standard ML bindings for Gelfgren *)

structure Gelfgren =
struct
    exception GelfgrenError of string

    (* Import C functions *)
    val gelfgren_last_error_message =
        _import "gelfgren_last_error_message" : unit -> MLton.Pointer.t;

    (* Bernstein Polynomial *)
    type bernstein_polynomial = MLton.Pointer.t

    val gelfgren_bernstein_create =
        _import "gelfgren_bernstein_create" : real array * int * real * real -> bernstein_polynomial;
    val gelfgren_bernstein_free =
        _import "gelfgren_bernstein_free" : bernstein_polynomial -> unit;
    val gelfgren_bernstein_evaluate =
        _import "gelfgren_bernstein_evaluate" : bernstein_polynomial * real * real ref -> int;
    val gelfgren_bernstein_derivative =
        _import "gelfgren_bernstein_derivative" : bernstein_polynomial -> bernstein_polynomial;
    val gelfgren_bernstein_integral =
        _import "gelfgren_bernstein_integral" : bernstein_polynomial -> bernstein_polynomial;
    val gelfgren_bernstein_degree =
        _import "gelfgren_bernstein_degree" : bernstein_polynomial * int ref -> int;
    val gelfgren_bernstein_add =
        _import "gelfgren_bernstein_add" : bernstein_polynomial * bernstein_polynomial -> bernstein_polynomial;
    val gelfgren_bernstein_scale =
        _import "gelfgren_bernstein_scale" : real * bernstein_polynomial -> bernstein_polynomial;

    (* Rational Function *)
    type rational_function = MLton.Pointer.t

    val gelfgren_rational_create =
        _import "gelfgren_rational_create" : bernstein_polynomial * bernstein_polynomial -> rational_function;
    val gelfgren_rational_free =
        _import "gelfgren_rational_free" : rational_function -> unit;
    val gelfgren_rational_evaluate =
        _import "gelfgren_rational_evaluate" : rational_function * real * real ref -> int;

    (* Padé Approximant *)
    type pade_approximant = MLton.Pointer.t

    val gelfgren_pade_from_series =
        _import "gelfgren_pade_from_series" : real array * int * int * int * real * real -> pade_approximant;
    val gelfgren_pade_free =
        _import "gelfgren_pade_free" : pade_approximant -> unit;
    val gelfgren_pade_evaluate =
        _import "gelfgren_pade_evaluate" : pade_approximant * real * real ref -> int;

    (* Helper functions *)
    fun getLastError () =
        let
            val ptr = gelfgren_last_error_message ()
        in
            if MLton.Pointer.isNull ptr then "Unknown error"
            else MLton.Pointer.getCString ptr
        end

    fun checkError (code, operation) =
        if code <> 0 then
            raise GelfgrenError (operation ^ ": " ^ getLastError ())
        else ()

    fun checkNull (ptr, operation) =
        if MLton.Pointer.isNull ptr then
            raise GelfgrenError (operation ^ ": " ^ getLastError ())
        else ptr

    (* High-level Bernstein Polynomial API *)
    structure BernsteinPolynomial =
    struct
        type t = bernstein_polynomial

        fun create (coeffs : real list, a : real, b : real) : t =
            let
                val arr = Array.fromList coeffs
                val degree = length coeffs - 1
                val ptr = gelfgren_bernstein_create (arr, degree, a, b)
            in
                checkNull (ptr, "create_bernstein")
            end

        fun free (poly : t) = gelfgren_bernstein_free poly

        fun evaluate (poly : t, x : real) : real =
            let
                val result = ref 0.0
                val code = gelfgren_bernstein_evaluate (poly, x, result)
                val _ = checkError (code, "evaluate")
            in
                !result
            end

        fun evaluateList (poly : t, xs : real list) : real list =
            map (fn x => evaluate (poly, x)) xs

        fun derivative (poly : t) : t =
            let
                val ptr = gelfgren_bernstein_derivative poly
            in
                checkNull (ptr, "derivative")
            end

        fun integral (poly : t) : t =
            let
                val ptr = gelfgren_bernstein_integral poly
            in
                checkNull (ptr, "integral")
            end

        fun degree (poly : t) : int =
            let
                val deg = ref 0
                val code = gelfgren_bernstein_degree (poly, deg)
                val _ = checkError (code, "degree")
            in
                !deg
            end

        fun add (p1 : t, p2 : t) : t =
            let
                val ptr = gelfgren_bernstein_add (p1, p2)
            in
                checkNull (ptr, "add")
            end

        fun scale (scalar : real, poly : t) : t =
            let
                val ptr = gelfgren_bernstein_scale (scalar, poly)
            in
                checkNull (ptr, "scale")
            end
    end

    (* High-level Rational Function API *)
    structure RationalFunction =
    struct
        type t = rational_function

        fun create (num : bernstein_polynomial, den : bernstein_polynomial) : t =
            let
                val ptr = gelfgren_rational_create (num, den)
            in
                checkNull (ptr, "create_rational")
            end

        fun free (rat : t) = gelfgren_rational_free rat

        fun evaluate (rat : t, x : real) : real =
            let
                val result = ref 0.0
                val code = gelfgren_rational_evaluate (rat, x, result)
                val _ = checkError (code, "evaluate_rational")
            in
                !result
            end

        fun evaluateList (rat : t, xs : real list) : real list =
            map (fn x => evaluate (rat, x)) xs
    end

    (* High-level Padé Approximant API *)
    structure PadeApproximant =
    struct
        type t = pade_approximant

        fun fromSeries (coeffs : real list, n : int, m : int, a : real, b : real) : t =
            let
                val arr = Array.fromList coeffs
                val len = length coeffs
                val ptr = gelfgren_pade_from_series (arr, len, n, m, a, b)
            in
                checkNull (ptr, "create_pade")
            end

        fun free (pade : t) = gelfgren_pade_free pade

        fun evaluate (pade : t, x : real) : real =
            let
                val result = ref 0.0
                val code = gelfgren_pade_evaluate (pade, x, result)
                val _ = checkError (code, "evaluate_pade")
            in
                !result
            end

        fun evaluateList (pade : t, xs : real list) : real list =
            map (fn x => evaluate (pade, x)) xs
    end
end
