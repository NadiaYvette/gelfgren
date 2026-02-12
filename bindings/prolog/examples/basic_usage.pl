:- use_module(library(gelfgren)).

main :-
    format('Gelfgren SWI-Prolog Bindings Example~n'),
    format('============================================================~n~n'),

    % Example 1: Bernstein Polynomial
    format('1. Bernstein Polynomial~n'),
    format('------------------------------------------------------------~n'),

    bernstein_create([1.0, 2.0, 6.0], 0.0, 1.0, Poly),
    bernstein_degree(Poly, Deg),
    format('Created polynomial with degree: ~d~n', [Deg]),

    bernstein_evaluate(Poly, 0.5, Y),
    format('f(0.5) = ~f~n', [Y]),

    % Vectorized evaluation
    format('~nVectorized evaluation:~n'),
    forall(member(X, [0.0, 0.25, 0.5, 0.75, 1.0]),
           (bernstein_evaluate(Poly, X, YVal),
            format('  f(~f) = ~f~n', [X, YVal]))),

    % Derivative
    format('~n'),
    bernstein_derivative(Poly, DPoly),
    bernstein_evaluate(DPoly, 0.5, DY),
    format('Derivative at x=0.5: ~f~n', [DY]),

    % Integral
    bernstein_integral(Poly, IPoly),
    bernstein_evaluate(IPoly, 1.0, Y1),
    bernstein_evaluate(IPoly, 0.0, Y0),
    Area is Y1 - Y0,
    format('Integral from 0 to 1: ~f~n', [Area]),

    % Clean up
    bernstein_free(IPoly),
    bernstein_free(DPoly),
    bernstein_free(Poly),

    format('~n'),

    % Example 2: Rational Function
    format('2. Rational Function~n'),
    format('------------------------------------------------------------~n'),

    bernstein_create([0.0, 1.0], 0.0, 1.0, Num),
    bernstein_create([1.0, 2.0], 0.0, 1.0, Den),
    rational_create(Num, Den, Rat),

    X = 0.5,
    rational_evaluate(Rat, X, YRat),
    Expected is X / (1.0 + X),
    format('R(0.5) = ~f~n', [YRat]),
    format('Expected: ~f~n', [Expected]),

    % Vectorized
    format('~nVectorized evaluation:~n'),
    forall(member(XRat, [0.0, 0.25, 0.5, 0.75, 1.0]),
           (rational_evaluate(Rat, XRat, YRatVal),
            Exp is XRat / (1.0 + XRat),
            format('  R(~f) = ~f (expected: ~f)~n', [XRat, YRatVal, Exp]))),

    % Clean up
    rational_free(Rat),
    bernstein_free(Den),
    bernstein_free(Num),

    format('~n'),
    format('All examples completed successfully!~n').

:- initialization(main, main).
