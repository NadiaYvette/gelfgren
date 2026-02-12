:- module(gelfgren, [
    bernstein_create/4,
    bernstein_free/1,
    bernstein_evaluate/3,
    bernstein_derivative/2,
    bernstein_integral/2,
    bernstein_degree/2,
    bernstein_add/3,
    bernstein_scale/3,
    rational_create/3,
    rational_free/1,
    rational_evaluate/3,
    pade_from_series/6,
    pade_free/1,
    pade_evaluate/3
]).

:- use_foreign_library(foreign(gelfgren)).

%! bernstein_create(+Coeffs:list(float), +A:float, +B:float, -Poly) is det.
%
%  Create a Bernstein polynomial from coefficients on interval [A, B].

bernstein_create(Coeffs, A, B, Poly) :-
    length(Coeffs, Len),
    Degree is Len - 1,
    bernstein_create_ffi(Coeffs, Degree, A, B, Poly).

:- foreign(bernstein_create_ffi, c,
           bernstein_create_ffi(+term, +integer, +float, +float, -address)).

%! bernstein_free(+Poly) is det.
%
%  Free a Bernstein polynomial.

:- foreign(bernstein_free, c, bernstein_free(+address)).

%! bernstein_evaluate(+Poly, +X:float, -Y:float) is det.
%
%  Evaluate Bernstein polynomial at point X.

bernstein_evaluate(Poly, X, Y) :-
    bernstein_evaluate_ffi(Poly, X, Y, Code),
    (Code =\= 0 -> throw(gelfgren_error('Evaluation failed')) ; true).

:- foreign(bernstein_evaluate_ffi, c,
           bernstein_evaluate_ffi(+address, +float, -float, [-integer])).

%! bernstein_derivative(+Poly, -DPoly) is det.
%
%  Compute the derivative of a Bernstein polynomial.

:- foreign(bernstein_derivative, c, bernstein_derivative(+address, -address)).

%! bernstein_integral(+Poly, -IPoly) is det.
%
%  Compute the antiderivative of a Bernstein polynomial.

:- foreign(bernstein_integral, c, bernstein_integral(+address, -address)).

%! bernstein_degree(+Poly, -Degree:integer) is det.
%
%  Get the degree of a Bernstein polynomial.

bernstein_degree(Poly, Degree) :-
    bernstein_degree_ffi(Poly, Degree, Code),
    (Code =\= 0 -> throw(gelfgren_error('Degree query failed')) ; true).

:- foreign(bernstein_degree_ffi, c,
           bernstein_degree_ffi(+address, -integer, [-integer])).

%! bernstein_add(+P1, +P2, -Sum) is det.
%
%  Add two Bernstein polynomials.

:- foreign(bernstein_add, c, bernstein_add(+address, +address, -address)).

%! bernstein_scale(+Scalar:float, +Poly, -Scaled) is det.
%
%  Scale a Bernstein polynomial by a scalar.

:- foreign(bernstein_scale, c, bernstein_scale(+float, +address, -address)).

%! rational_create(+Num, +Den, -Rat) is det.
%
%  Create a rational function from numerator and denominator.

:- foreign(rational_create, c, rational_create(+address, +address, -address)).

%! rational_free(+Rat) is det.
%
%  Free a rational function.

:- foreign(rational_free, c, rational_free(+address)).

%! rational_evaluate(+Rat, +X:float, -Y:float) is det.
%
%  Evaluate a rational function at point X.

rational_evaluate(Rat, X, Y) :-
    rational_evaluate_ffi(Rat, X, Y, Code),
    (Code =\= 0 -> throw(gelfgren_error('Rational evaluation failed')) ; true).

:- foreign(rational_evaluate_ffi, c,
           rational_evaluate_ffi(+address, +float, -float, [-integer])).

%! pade_from_series(+Coeffs:list(float), +N:int, +M:int, +A:float, +B:float, -Pade) is det.
%
%  Create a Padé approximant [N/M] from power series coefficients.

pade_from_series(Coeffs, N, M, A, B, Pade) :-
    length(Coeffs, Len),
    pade_from_series_ffi(Coeffs, Len, N, M, A, B, Pade).

:- foreign(pade_from_series_ffi, c,
           pade_from_series_ffi(+term, +integer, +integer, +integer, +float, +float, -address)).

%! pade_free(+Pade) is det.
%
%  Free a Padé approximant.

:- foreign(pade_free, c, pade_free(+address)).

%! pade_evaluate(+Pade, +X:float, -Y:float) is det.
%
%  Evaluate a Padé approximant at point X.

pade_evaluate(Pade, X, Y) :-
    pade_evaluate_ffi(Pade, X, Y, Code),
    (Code =\= 0 -> throw(gelfgren_error('Padé evaluation failed')) ; true).

:- foreign(pade_evaluate_ffi, c,
           pade_evaluate_ffi(+address, +float, -float, [-integer])).
