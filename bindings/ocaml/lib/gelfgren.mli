(** OCaml bindings for Gelfgren piecewise rational interpolation library *)

(** Exception raised for Gelfgren errors *)
exception Gelfgren_error of string

(** Bernstein polynomial operations *)
module BernsteinPolynomial : sig
  (** Type representing a Bernstein polynomial *)
  type t

  (** [create coeffs a b] creates a Bernstein polynomial with given
      coefficients on interval \[a, b\].
      @raise Gelfgren_error if coefficients array is empty or interval is invalid *)
  val create : float array -> float -> float -> t

  (** [free poly] frees the memory associated with a polynomial.
      Safe to call multiple times. *)
  val free : t -> unit

  (** [with_poly coeffs a b f] creates a polynomial, applies function [f] to it,
      and ensures cleanup even if [f] raises an exception. *)
  val with_poly : float array -> float -> float -> (t -> 'a) -> 'a

  (** [evaluate poly x] evaluates the polynomial at point [x].
      @raise Gelfgren_error if polynomial has been freed or evaluation fails *)
  val evaluate : t -> float -> float

  (** [evaluate_array poly xs] evaluates the polynomial at multiple points.
      Equivalent to [Array.map (evaluate poly) xs] *)
  val evaluate_array : t -> float array -> float array

  (** [derivative poly] computes the derivative of the polynomial.
      Returns a new polynomial that must be freed separately.
      @raise Gelfgren_error on failure *)
  val derivative : t -> t

  (** [integral poly] computes the antiderivative of the polynomial.
      Returns a new polynomial that must be freed separately.
      @raise Gelfgren_error on failure *)
  val integral : t -> t

  (** [degree poly] returns the degree of the polynomial.
      @raise Gelfgren_error if polynomial has been freed *)
  val degree : t -> int

  (** [add p1 p2] adds two polynomials.
      Returns a new polynomial that must be freed separately.
      @raise Gelfgren_error if either polynomial has been freed *)
  val add : t -> t -> t

  (** [subtract p1 p2] subtracts [p2] from [p1].
      Returns a new polynomial that must be freed separately. *)
  val subtract : t -> t -> t

  (** [multiply p1 p2] multiplies two polynomials.
      Returns a new polynomial that must be freed separately. *)
  val multiply : t -> t -> t

  (** [scale scalar poly] multiplies polynomial by a scalar.
      Returns a new polynomial that must be freed separately. *)
  val scale : float -> t -> t

  (** Infix operator for addition *)
  val ( + ) : t -> t -> t

  (** Infix operator for subtraction *)
  val ( - ) : t -> t -> t

  (** Infix operator for multiplication *)
  val ( * ) : t -> t -> t

  (** Infix operator for scalar multiplication *)
  val ( *. ) : float -> t -> t
end

(** Rational function operations *)
module RationalFunction : sig
  (** Type representing a rational function P(x)/Q(x) *)
  type t

  (** [create num den] creates a rational function from numerator and
      denominator polynomials.
      @raise Gelfgren_error if either polynomial has been freed *)
  val create : BernsteinPolynomial.t -> BernsteinPolynomial.t -> t

  (** [free rat] frees the memory associated with a rational function *)
  val free : t -> unit

  (** [with_rational num den f] creates a rational function, applies [f],
      and ensures cleanup *)
  val with_rational : BernsteinPolynomial.t -> BernsteinPolynomial.t ->
    (t -> 'a) -> 'a

  (** [evaluate rat x] evaluates the rational function at point [x].
      @raise Gelfgren_error if denominator is zero (pole) *)
  val evaluate : t -> float -> float

  (** [evaluate_array rat xs] evaluates at multiple points *)
  val evaluate_array : t -> float array -> float array

  (** [derivative rat] computes the derivative using quotient rule.
      Returns a new rational function that must be freed separately. *)
  val derivative : t -> t
end

(** Padé approximant operations *)
module PadeApproximant : sig
  (** Type representing a Padé approximant \[n/m\] *)
  type t

  (** [from_series coeffs n m a b] creates a Padé approximant of order \[n/m\]
      from power series coefficients on interval \[a, b\].
      @param coeffs Power series coefficients \[c_0, c_1, ..., c_k\]
      @param n Numerator degree
      @param m Denominator degree
      @raise Gelfgren_error if construction fails (e.g., singular system) *)
  val from_series : float array -> int -> int -> float -> float -> t

  (** [free pade] frees the memory associated with a Padé approximant *)
  val free : t -> unit

  (** [with_pade coeffs n m a b f] creates a Padé approximant, applies [f],
      and ensures cleanup *)
  val with_pade : float array -> int -> int -> float -> float ->
    (t -> 'a) -> 'a

  (** [evaluate pade x] evaluates the Padé approximant at point [x].
      @raise Gelfgren_error if denominator is zero (pole) *)
  val evaluate : t -> float -> float

  (** [evaluate_array pade xs] evaluates at multiple points *)
  val evaluate_array : t -> float array -> float array
end
