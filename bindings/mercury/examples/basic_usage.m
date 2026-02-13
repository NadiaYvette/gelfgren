%---------------------------------------------------------------------------%
% vim: ft=mercury ts=4 sw=4 et
%---------------------------------------------------------------------------%
% File: basic_usage.m
% Main author: Gelfgren Contributors
%
% Example usage of Gelfgren Mercury bindings.
%
%---------------------------------------------------------------------------%

:- module basic_usage.
:- interface.

:- import_module io.

:- pred main(io::di, io::uo) is det.

%---------------------------------------------------------------------------%
%---------------------------------------------------------------------------%

:- implementation.

:- import_module float.
:- import_module gelfgren.
:- import_module list.
:- import_module math.
:- import_module maybe.
:- import_module string.

%---------------------------------------------------------------------------%

main(!IO) :-
    io.write_string("Gelfgren Mercury Bindings Example\n", !IO),
    io.write_string(string.duplicate_char('=', 60) ++ "\n\n", !IO),

    example1(!IO),
    example2(!IO),
    example3(!IO),
    example4(!IO),

    io.nl(!IO),
    io.write_string("All examples completed successfully!\n", !IO).

%---------------------------------------------------------------------------%
% Example 1: Bernstein Polynomial
%---------------------------------------------------------------------------%

:- pred example1(io::di, io::uo) is det.

example1(!IO) :-
    io.write_string("1. Bernstein Polynomial\n", !IO),
    io.write_string(string.duplicate_char('-', 60) ++ "\n", !IO),

    % Create a quadratic polynomial: 2x^2 + 3x + 1 on [0, 1]
    % Bernstein coefficients: [1, 2, 6]
    Coeffs = [1.0, 2.0, 6.0],
    create_bernstein(Coeffs, 0.0, 1.0, PolyResult, !IO),

    (
        PolyResult = ok(Poly),

        % Get degree
        degree_bernstein(Poly, DegResult, !IO),
        (
            DegResult = ok(Deg),
            io.format("Created polynomial with degree: %d\n", [i(Deg)], !IO)
        ;
            DegResult = error(Err),
            io.format("Error getting degree: %s\n",
                [s(error_to_string(Err))], !IO)
        ),

        % Evaluate at a point
        evaluate_bernstein(Poly, 0.5, EvalResult, !IO),
        (
            EvalResult = ok(Y),
            io.format("f(0.5) = %.6f\n", [f(Y)], !IO)
        ;
            EvalResult = error(Err),
            io.format("Error evaluating: %s\n",
                [s(error_to_string(Err))], !IO)
        ),

        % Vectorized evaluation
        io.nl(!IO),
        io.write_string("Vectorized evaluation:\n", !IO),
        XValues = [0.0, 0.25, 0.5, 0.75, 1.0],
        list.foldl(evaluate_and_print(Poly), XValues, !IO),

        % Derivative
        io.nl(!IO),
        derivative_bernstein(Poly, DPolyResult, !IO),
        (
            DPolyResult = ok(DPoly),
            evaluate_bernstein(DPoly, 0.5, DEvalResult, !IO),
            (
                DEvalResult = ok(DY),
                io.format("Derivative at x=0.5: %.6f\n", [f(DY)], !IO)
            ;
                DEvalResult = error(Err),
                io.format("Error: %s\n", [s(error_to_string(Err))], !IO)
            ),
            free_bernstein(DPoly, !IO)
        ;
            DPolyResult = error(Err),
            io.format("Error computing derivative: %s\n",
                [s(error_to_string(Err))], !IO)
        ),

        % Integral
        integral_bernstein(Poly, IPolyResult, !IO),
        (
            IPolyResult = ok(IPoly),
            evaluate_bernstein(IPoly, 1.0, IEval1Result, !IO),
            evaluate_bernstein(IPoly, 0.0, IEval0Result, !IO),
            (
                IEval1Result = ok(Y1),
                IEval0Result = ok(Y0)
            ->
                Area = Y1 - Y0,
                io.format("Integral from 0 to 1: %.6f\n", [f(Area)], !IO)
            ;
                io.write_string("Error evaluating integral\n", !IO)
            ),
            free_bernstein(IPoly, !IO)
        ;
            IPolyResult = error(Err),
            io.format("Error computing integral: %s\n",
                [s(error_to_string(Err))], !IO)
        ),

        free_bernstein(Poly, !IO)
    ;
        PolyResult = error(Err),
        io.format("Error creating polynomial: %s\n",
            [s(error_to_string(Err))], !IO)
    ),

    io.nl(!IO).

:- pred evaluate_and_print(bernstein_polynomial::in, float::in,
    io::di, io::uo) is det.

evaluate_and_print(Poly, X, !IO) :-
    evaluate_bernstein(Poly, X, Result, !IO),
    (
        Result = ok(Y),
        io.format("  f(%.2f) = %.6f\n", [f(X), f(Y)], !IO)
    ;
        Result = error(Err),
        io.format("  Error at %.2f: %s\n",
            [f(X), s(error_to_string(Err))], !IO)
    ).

%---------------------------------------------------------------------------%
% Example 2: Polynomial Arithmetic
%---------------------------------------------------------------------------%

:- pred example2(io::di, io::uo) is det.

example2(!IO) :-
    io.write_string("2. Polynomial Arithmetic\n", !IO),
    io.write_string(string.duplicate_char('-', 60) ++ "\n", !IO),

    create_bernstein([1.0, 2.0, 3.0], 0.0, 1.0, P1Result, !IO),
    create_bernstein([2.0, 1.0, 1.0], 0.0, 1.0, P2Result, !IO),

    (
        P1Result = ok(P1),
        P2Result = ok(P2)
    ->
        evaluate_bernstein(P1, 0.5, Eval1, !IO),
        evaluate_bernstein(P2, 0.5, Eval2, !IO),

        (
            Eval1 = ok(Y1),
            Eval2 = ok(Y2)
        ->
            io.format("p1(0.5) = %.6f\n", [f(Y1)], !IO),
            io.format("p2(0.5) = %.6f\n", [f(Y2)], !IO)
        ;
            io.write_string("Error evaluating polynomials\n", !IO)
        ),

        % Addition
        add_bernstein(P1, P2, SumResult, !IO),
        (
            SumResult = ok(Sum),
            evaluate_bernstein(Sum, 0.5, SumEval, !IO),
            (
                SumEval = ok(YSum),
                io.format("p1 + p2 at x=0.5: %.6f\n", [f(YSum)], !IO)
            ;
                SumEval = error(_),
                io.write_string("Error evaluating sum\n", !IO)
            ),
            free_bernstein(Sum, !IO)
        ;
            SumResult = error(Err),
            io.format("Error adding: %s\n", [s(error_to_string(Err))], !IO)
        ),

        % Scaling
        scale_bernstein(2.0, P1, ScaledResult, !IO),
        (
            ScaledResult = ok(Scaled),
            evaluate_bernstein(Scaled, 0.5, ScaledEval, !IO),
            (
                ScaledEval = ok(YScaled),
                io.format("2 * p1 at x=0.5: %.6f\n", [f(YScaled)], !IO)
            ;
                ScaledEval = error(_),
                io.write_string("Error evaluating scaled\n", !IO)
            ),
            free_bernstein(Scaled, !IO)
        ;
            ScaledResult = error(Err),
            io.format("Error scaling: %s\n", [s(error_to_string(Err))], !IO)
        ),

        free_bernstein(P1, !IO),
        free_bernstein(P2, !IO)
    ;
        io.write_string("Error creating polynomials\n", !IO)
    ),

    io.nl(!IO).

%---------------------------------------------------------------------------%
% Example 3: Rational Function
%---------------------------------------------------------------------------%

:- pred example3(io::di, io::uo) is det.

example3(!IO) :-
    io.write_string("3. Rational Function\n", !IO),
    io.write_string(string.duplicate_char('-', 60) ++ "\n", !IO),

    % Create rational function: x / (1 + x) on [0, 1]
    create_bernstein([0.0, 1.0], 0.0, 1.0, NumResult, !IO),
    create_bernstein([1.0, 2.0], 0.0, 1.0, DenResult, !IO),

    (
        NumResult = ok(Num),
        DenResult = ok(Den)
    ->
        create_rational(Num, Den, RatResult, !IO),
        (
            RatResult = ok(Rat),

            % Evaluate
            X = 0.5,
            evaluate_rational(Rat, X, EvalResult, !IO),
            (
                EvalResult = ok(Y),
                Expected = X / (1.0 + X),
                io.format("R(%.1f) = %.6f\n", [f(X), f(Y)], !IO),
                io.format("Expected: %.6f\n", [f(Expected)], !IO)
            ;
                EvalResult = error(Err),
                io.format("Error: %s\n", [s(error_to_string(Err))], !IO)
            ),

            % Vectorized evaluation
            io.nl(!IO),
            io.write_string("Vectorized evaluation:\n", !IO),
            XValues = [0.0, 0.25, 0.5, 0.75, 1.0],
            list.foldl(evaluate_rational_and_print(Rat), XValues, !IO),

            free_rational(Rat, !IO)
        ;
            RatResult = error(Err),
            io.format("Error creating rational: %s\n",
                [s(error_to_string(Err))], !IO)
        ),

        free_bernstein(Num, !IO),
        free_bernstein(Den, !IO)
    ;
        io.write_string("Error creating numerator/denominator\n", !IO)
    ),

    io.nl(!IO).

:- pred evaluate_rational_and_print(rational_function::in, float::in,
    io::di, io::uo) is det.

evaluate_rational_and_print(Rat, X, !IO) :-
    evaluate_rational(Rat, X, Result, !IO),
    (
        Result = ok(Y),
        Expected = X / (1.0 + X),
        io.format("  R(%.2f) = %.6f (expected: %.6f)\n",
            [f(X), f(Y), f(Expected)], !IO)
    ;
        Result = error(Err),
        io.format("  Error at %.2f: %s\n",
            [f(X), s(error_to_string(Err))], !IO)
    ).

%---------------------------------------------------------------------------%
% Example 4: Padé Approximant
%---------------------------------------------------------------------------%

:- pred example4(io::di, io::uo) is det.

example4(!IO) :-
    io.write_string("4. Padé Approximant\n", !IO),
    io.write_string(string.duplicate_char('-', 60) ++ "\n", !IO),

    % Approximate exp(x) near x=0 with [2/2] Padé approximant
    % Power series: 1 + x + x^2/2 + x^3/6 + x^4/24
    Coeffs = [1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0],
    create_pade_from_series(Coeffs, 2, 2, -1.0, 1.0, PadeResult, !IO),

    (
        PadeResult = ok(Pade),

        io.write_string("Approximation of exp(x):\n", !IO),
        XValues = [-0.5, -0.25, 0.0, 0.25, 0.5],
        list.foldl(evaluate_pade_and_print(Pade), XValues, !IO),

        free_pade(Pade, !IO)
    ;
        PadeResult = error(Err),
        io.format("Error creating Padé: %s\n",
            [s(error_to_string(Err))], !IO)
    ),

    io.nl(!IO).

:- pred evaluate_pade_and_print(pade_approximant::in, float::in,
    io::di, io::uo) is det.

evaluate_pade_and_print(Pade, X, !IO) :-
    evaluate_pade(Pade, X, Result, !IO),
    (
        Result = ok(Y),
        Exact = math.exp(X),
        Error = float.abs(Y - Exact),
        io.format("  x=%.2f: Padé=%.6f, exp(x)=%.6f, error=%.8f\n",
            [f(X), f(Y), f(Exact), f(Error)], !IO)
    ;
        Result = error(Err),
        io.format("  Error at x=%.2f: %s\n",
            [f(X), s(error_to_string(Err))], !IO)
    ).

%---------------------------------------------------------------------------%
:- end_module basic_usage.
%---------------------------------------------------------------------------%
