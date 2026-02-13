%---------------------------------------------------------------------------%
% vim: ft=mercury ts=4 sw=4 et
%---------------------------------------------------------------------------%
% File: gelfgren.m
% Main authors: Gelfgren Contributors
% Stability: low
%
% Mercury bindings for the Gelfgren numerical computing library.
%
% This module provides Mercury bindings to the Gelfgren library, which
% implements piecewise rational interpolation methods based on Jan Gelfgren's
% 1975 research.
%
%---------------------------------------------------------------------------%

:- module gelfgren.
:- interface.

:- import_module io.
:- import_module list.
:- import_module maybe.
:- import_module string.

%---------------------------------------------------------------------------%
% Types
%---------------------------------------------------------------------------%

    % Opaque type for Bernstein polynomial
:- type bernstein_polynomial.

    % Opaque type for rational function
:- type rational_function.

    % Opaque type for Padé approximant
:- type pade_approximant.

    % Error type
:- type gelfgren_error
    --->    null_pointer_error(string)
    ;       invalid_interval_error(string)
    ;       empty_coefficients_error(string)
    ;       singular_matrix_error(string)
    ;       division_by_zero_error(string)
    ;       pole_error(string)
    ;       unknown_error(string).

    % Result type for operations that may fail
:- type gelfgren_result(T) == maybe_error(T, gelfgren_error).

%---------------------------------------------------------------------------%
% Bernstein polynomial operations
%---------------------------------------------------------------------------%

    % Create a Bernstein polynomial from coefficients on interval [a, b]
    %
:- pred create_bernstein(list(float)::in, float::in, float::in,
    gelfgren_result(bernstein_polynomial)::out, io::di, io::uo) is det.

    % Free a Bernstein polynomial
    %
:- pred free_bernstein(bernstein_polynomial::in, io::di, io::uo) is det.

    % Evaluate a Bernstein polynomial at a point
    %
:- pred evaluate_bernstein(bernstein_polynomial::in, float::in,
    gelfgren_result(float)::out, io::di, io::uo) is det.

    % Compute the derivative of a Bernstein polynomial
    %
:- pred derivative_bernstein(bernstein_polynomial::in,
    gelfgren_result(bernstein_polynomial)::out, io::di, io::uo) is det.

    % Compute the integral of a Bernstein polynomial
    %
:- pred integral_bernstein(bernstein_polynomial::in,
    gelfgren_result(bernstein_polynomial)::out, io::di, io::uo) is det.

    % Get the degree of a Bernstein polynomial
    %
:- pred degree_bernstein(bernstein_polynomial::in,
    gelfgren_result(int)::out, io::di, io::uo) is det.

    % Add two Bernstein polynomials
    %
:- pred add_bernstein(bernstein_polynomial::in, bernstein_polynomial::in,
    gelfgren_result(bernstein_polynomial)::out, io::di, io::uo) is det.

    % Subtract two Bernstein polynomials
    %
:- pred subtract_bernstein(bernstein_polynomial::in, bernstein_polynomial::in,
    gelfgren_result(bernstein_polynomial)::out, io::di, io::uo) is det.

    % Multiply two Bernstein polynomials
    %
:- pred multiply_bernstein(bernstein_polynomial::in, bernstein_polynomial::in,
    gelfgren_result(bernstein_polynomial)::out, io::di, io::uo) is det.

    % Scale a Bernstein polynomial by a scalar
    %
:- pred scale_bernstein(float::in, bernstein_polynomial::in,
    gelfgren_result(bernstein_polynomial)::out, io::di, io::uo) is det.

%---------------------------------------------------------------------------%
% Rational function operations
%---------------------------------------------------------------------------%

    % Create a rational function from numerator and denominator polynomials
    %
:- pred create_rational(bernstein_polynomial::in, bernstein_polynomial::in,
    gelfgren_result(rational_function)::out, io::di, io::uo) is det.

    % Free a rational function
    %
:- pred free_rational(rational_function::in, io::di, io::uo) is det.

    % Evaluate a rational function at a point
    %
:- pred evaluate_rational(rational_function::in, float::in,
    gelfgren_result(float)::out, io::di, io::uo) is det.

    % Compute the derivative of a rational function
    %
:- pred derivative_rational(rational_function::in,
    gelfgren_result(rational_function)::out, io::di, io::uo) is det.

%---------------------------------------------------------------------------%
% Padé approximant operations
%---------------------------------------------------------------------------%

    % Create a Padé approximant from power series coefficients
    %
:- pred create_pade_from_series(list(float)::in, int::in, int::in,
    float::in, float::in, gelfgren_result(pade_approximant)::out,
    io::di, io::uo) is det.

    % Free a Padé approximant
    %
:- pred free_pade(pade_approximant::in, io::di, io::uo) is det.

    % Evaluate a Padé approximant at a point
    %
:- pred evaluate_pade(pade_approximant::in, float::in,
    gelfgren_result(float)::out, io::di, io::uo) is det.

%---------------------------------------------------------------------------%
% Utility predicates
%---------------------------------------------------------------------------%

    % Convert error to string for display
    %
:- func error_to_string(gelfgren_error) = string.

%---------------------------------------------------------------------------%
%---------------------------------------------------------------------------%

:- implementation.

:- import_module bool.
:- import_module float.
:- import_module int.
:- import_module require.

%---------------------------------------------------------------------------%
% C foreign type declarations
%---------------------------------------------------------------------------%

:- pragma foreign_type("C", bernstein_polynomial, "void *").
:- pragma foreign_type("C", rational_function, "void *").
:- pragma foreign_type("C", pade_approximant, "void *").

%---------------------------------------------------------------------------%
% Error handling
%---------------------------------------------------------------------------%

:- pred get_last_error(string::out, io::di, io::uo) is det.

:- pragma foreign_proc("C",
    get_last_error(Msg::out, _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    const char *c_msg = gelfgren_last_error_message();
    if (c_msg != NULL) {
        MR_make_aligned_string_copy(Msg, c_msg);
    } else {
        MR_make_aligned_string_copy(Msg, \"Unknown error\");
    }
").

error_to_string(null_pointer_error(Msg)) = "Null pointer: " ++ Msg.
error_to_string(invalid_interval_error(Msg)) = "Invalid interval: " ++ Msg.
error_to_string(empty_coefficients_error(Msg)) = "Empty coefficients: " ++ Msg.
error_to_string(singular_matrix_error(Msg)) = "Singular matrix: " ++ Msg.
error_to_string(division_by_zero_error(Msg)) = "Division by zero: " ++ Msg.
error_to_string(pole_error(Msg)) = "Pole detected: " ++ Msg.
error_to_string(unknown_error(Msg)) = "Unknown error: " ++ Msg.

%---------------------------------------------------------------------------%
% Bernstein polynomial operations
%---------------------------------------------------------------------------%

create_bernstein(Coeffs, A, B, Result, !IO) :-
    list.length(Coeffs, Len),
    ( if Len = 0 then
        Result = error(empty_coefficients_error("Coefficient list is empty"))
    else
        Degree = Len - 1,
        c_bernstein_create(Coeffs, Degree, A, B, Ptr, !IO),
        ( if is_null_ptr(Ptr) then
            get_last_error(Msg, !IO),
            Result = error(null_pointer_error(Msg))
        else
            Result = ok(Ptr)
        )
    ).

:- pred c_bernstein_create(list(float)::in, int::in, float::in, float::in,
    bernstein_polynomial::out, io::di, io::uo) is det.

:- pragma foreign_proc("C",
    c_bernstein_create(Coeffs::in, Degree::in, A::in, B::in, Ptr::out,
        _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    int n = 0;
    double *coeffs_array;
    MR_Word list = Coeffs;

    // Count elements
    while (!MR_list_is_empty(list)) {
        n++;
        list = MR_list_tail(list);
    }

    // Allocate array
    coeffs_array = MR_GC_NEW_ARRAY(double, n);

    // Copy coefficients
    list = Coeffs;
    for (int i = 0; i < n; i++) {
        coeffs_array[i] = MR_word_to_float(MR_list_head(list));
        list = MR_list_tail(list);
    }

    Ptr = (MR_Word)gelfgren_bernstein_create(coeffs_array, (size_t)Degree, A, B);
").

free_bernstein(Ptr, !IO) :-
    c_bernstein_free(Ptr, !IO).

:- pragma foreign_proc("C",
    c_bernstein_free(Ptr::in, _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    if (Ptr != NULL) {
        gelfgren_bernstein_free((void *)Ptr);
    }
").

evaluate_bernstein(Ptr, X, Result, !IO) :-
    c_bernstein_evaluate(Ptr, X, Success, Y, !IO),
    ( if Success = yes then
        Result = ok(Y)
    else
        get_last_error(Msg, !IO),
        Result = error(unknown_error(Msg))
    ).

:- pred c_bernstein_evaluate(bernstein_polynomial::in, float::in, bool::out,
    float::out, io::di, io::uo) is det.

:- pragma foreign_proc("C",
    c_bernstein_evaluate(Ptr::in, X::in, Success::out, Result::out,
        _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    double result;
    int code = gelfgren_bernstein_evaluate((void *)Ptr, X, &result);
    Success = (code == 0) ? MR_YES : MR_NO;
    Result = result;
").

derivative_bernstein(Ptr, Result, !IO) :-
    c_bernstein_derivative(Ptr, DPtr, !IO),
    ( if is_null_ptr(DPtr) then
        get_last_error(Msg, !IO),
        Result = error(null_pointer_error(Msg))
    else
        Result = ok(DPtr)
    ).

:- pred c_bernstein_derivative(bernstein_polynomial::in,
    bernstein_polynomial::out, io::di, io::uo) is det.

:- pragma foreign_proc("C",
    c_bernstein_derivative(Ptr::in, DPtr::out, _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    DPtr = (MR_Word)gelfgren_bernstein_derivative((void *)Ptr);
").

integral_bernstein(Ptr, Result, !IO) :-
    c_bernstein_integral(Ptr, IPtr, !IO),
    ( if is_null_ptr(IPtr) then
        get_last_error(Msg, !IO),
        Result = error(null_pointer_error(Msg))
    else
        Result = ok(IPtr)
    ).

:- pred c_bernstein_integral(bernstein_polynomial::in,
    bernstein_polynomial::out, io::di, io::uo) is det.

:- pragma foreign_proc("C",
    c_bernstein_integral(Ptr::in, IPtr::out, _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    IPtr = (MR_Word)gelfgren_bernstein_integral((void *)Ptr);
").

degree_bernstein(Ptr, Result, !IO) :-
    c_bernstein_degree(Ptr, Success, Deg, !IO),
    ( if Success = yes then
        Result = ok(Deg)
    else
        get_last_error(Msg, !IO),
        Result = error(unknown_error(Msg))
    ).

:- pred c_bernstein_degree(bernstein_polynomial::in, bool::out, int::out,
    io::di, io::uo) is det.

:- pragma foreign_proc("C",
    c_bernstein_degree(Ptr::in, Success::out, Degree::out,
        _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    size_t deg;
    int code = gelfgren_bernstein_degree((void *)Ptr, &deg);
    Success = (code == 0) ? MR_YES : MR_NO;
    Degree = (MR_Integer)deg;
").

add_bernstein(Ptr1, Ptr2, Result, !IO) :-
    c_bernstein_add(Ptr1, Ptr2, SPtr, !IO),
    ( if is_null_ptr(SPtr) then
        get_last_error(Msg, !IO),
        Result = error(null_pointer_error(Msg))
    else
        Result = ok(SPtr)
    ).

:- pred c_bernstein_add(bernstein_polynomial::in, bernstein_polynomial::in,
    bernstein_polynomial::out, io::di, io::uo) is det.

:- pragma foreign_proc("C",
    c_bernstein_add(Ptr1::in, Ptr2::in, SPtr::out, _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    SPtr = (MR_Word)gelfgren_bernstein_add((void *)Ptr1, (void *)Ptr2);
").

subtract_bernstein(Ptr1, Ptr2, Result, !IO) :-
    c_bernstein_subtract(Ptr1, Ptr2, DPtr, !IO),
    ( if is_null_ptr(DPtr) then
        get_last_error(Msg, !IO),
        Result = error(null_pointer_error(Msg))
    else
        Result = ok(DPtr)
    ).

:- pred c_bernstein_subtract(bernstein_polynomial::in, bernstein_polynomial::in,
    bernstein_polynomial::out, io::di, io::uo) is det.

:- pragma foreign_proc("C",
    c_bernstein_subtract(Ptr1::in, Ptr2::in, DPtr::out, _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    DPtr = (MR_Word)gelfgren_bernstein_subtract((void *)Ptr1, (void *)Ptr2);
").

multiply_bernstein(Ptr1, Ptr2, Result, !IO) :-
    c_bernstein_multiply(Ptr1, Ptr2, PPtr, !IO),
    ( if is_null_ptr(PPtr) then
        get_last_error(Msg, !IO),
        Result = error(null_pointer_error(Msg))
    else
        Result = ok(PPtr)
    ).

:- pred c_bernstein_multiply(bernstein_polynomial::in, bernstein_polynomial::in,
    bernstein_polynomial::out, io::di, io::uo) is det.

:- pragma foreign_proc("C",
    c_bernstein_multiply(Ptr1::in, Ptr2::in, PPtr::out, _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    PPtr = (MR_Word)gelfgren_bernstein_multiply((void *)Ptr1, (void *)Ptr2);
").

scale_bernstein(Scalar, Ptr, Result, !IO) :-
    c_bernstein_scale(Ptr, Scalar, SPtr, !IO),
    ( if is_null_ptr(SPtr) then
        get_last_error(Msg, !IO),
        Result = error(null_pointer_error(Msg))
    else
        Result = ok(SPtr)
    ).

:- pred c_bernstein_scale(bernstein_polynomial::in, float::in,
    bernstein_polynomial::out, io::di, io::uo) is det.

:- pragma foreign_proc("C",
    c_bernstein_scale(Ptr::in, Scalar::in, SPtr::out, _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    SPtr = (MR_Word)gelfgren_bernstein_scale((void *)Ptr, Scalar);
").

%---------------------------------------------------------------------------%
% Rational function operations
%---------------------------------------------------------------------------%

create_rational(NumPtr, DenPtr, Result, !IO) :-
    c_rational_create(NumPtr, DenPtr, RPtr, !IO),
    ( if is_null_ptr(RPtr) then
        get_last_error(Msg, !IO),
        Result = error(null_pointer_error(Msg))
    else
        Result = ok(RPtr)
    ).

:- pred c_rational_create(bernstein_polynomial::in, bernstein_polynomial::in,
    rational_function::out, io::di, io::uo) is det.

:- pragma foreign_proc("C",
    c_rational_create(NumPtr::in, DenPtr::in, RPtr::out, _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    RPtr = (MR_Word)gelfgren_rational_create((void *)NumPtr, (void *)DenPtr);
").

free_rational(Ptr, !IO) :-
    c_rational_free(Ptr, !IO).

:- pragma foreign_proc("C",
    c_rational_free(Ptr::in, _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    if (Ptr != NULL) {
        gelfgren_rational_free((void *)Ptr);
    }
").

evaluate_rational(Ptr, X, Result, !IO) :-
    c_rational_evaluate(Ptr, X, Success, Y, !IO),
    ( if Success = yes then
        Result = ok(Y)
    else
        get_last_error(Msg, !IO),
        Result = error(pole_error(Msg))
    ).

:- pred c_rational_evaluate(rational_function::in, float::in, bool::out,
    float::out, io::di, io::uo) is det.

:- pragma foreign_proc("C",
    c_rational_evaluate(Ptr::in, X::in, Success::out, Result::out,
        _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    double result;
    int code = gelfgren_rational_evaluate((void *)Ptr, X, &result);
    Success = (code == 0) ? MR_YES : MR_NO;
    Result = result;
").

derivative_rational(Ptr, Result, !IO) :-
    c_rational_derivative(Ptr, DPtr, !IO),
    ( if is_null_ptr(DPtr) then
        get_last_error(Msg, !IO),
        Result = error(null_pointer_error(Msg))
    else
        Result = ok(DPtr)
    ).

:- pred c_rational_derivative(rational_function::in, rational_function::out,
    io::di, io::uo) is det.

:- pragma foreign_proc("C",
    c_rational_derivative(Ptr::in, DPtr::out, _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    DPtr = (MR_Word)gelfgren_rational_derivative((void *)Ptr);
").

%---------------------------------------------------------------------------%
% Padé approximant operations
%---------------------------------------------------------------------------%

create_pade_from_series(Coeffs, N, M, A, B, Result, !IO) :-
    list.length(Coeffs, Len),
    ( if Len = 0 then
        Result = error(empty_coefficients_error("Coefficient list is empty"))
    else
        c_pade_from_series(Coeffs, Len, N, M, A, B, PPtr, !IO),
        ( if is_null_ptr(PPtr) then
            get_last_error(Msg, !IO),
            Result = error(null_pointer_error(Msg))
        else
            Result = ok(PPtr)
        )
    ).

:- pred c_pade_from_series(list(float)::in, int::in, int::in, int::in,
    float::in, float::in, pade_approximant::out, io::di, io::uo) is det.

:- pragma foreign_proc("C",
    c_pade_from_series(Coeffs::in, NCoeffs::in, N::in, M::in, A::in, B::in,
        PPtr::out, _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    int n = 0;
    double *coeffs_array;
    MR_Word list = Coeffs;

    // Count elements
    while (!MR_list_is_empty(list)) {
        n++;
        list = MR_list_tail(list);
    }

    // Allocate array
    coeffs_array = MR_GC_NEW_ARRAY(double, n);

    // Copy coefficients
    list = Coeffs;
    for (int i = 0; i < n; i++) {
        coeffs_array[i] = MR_word_to_float(MR_list_head(list));
        list = MR_list_tail(list);
    }

    PPtr = (MR_Word)gelfgren_pade_from_power_series(
        coeffs_array, (size_t)NCoeffs, (size_t)N, (size_t)M, A, B);
").

free_pade(Ptr, !IO) :-
    c_pade_free(Ptr, !IO).

:- pragma foreign_proc("C",
    c_pade_free(Ptr::in, _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    if (Ptr != NULL) {
        gelfgren_pade_free((void *)Ptr);
    }
").

evaluate_pade(Ptr, X, Result, !IO) :-
    c_pade_evaluate(Ptr, X, Success, Y, !IO),
    ( if Success = yes then
        Result = ok(Y)
    else
        get_last_error(Msg, !IO),
        Result = error(pole_error(Msg))
    ).

:- pred c_pade_evaluate(pade_approximant::in, float::in, bool::out,
    float::out, io::di, io::uo) is det.

:- pragma foreign_proc("C",
    c_pade_evaluate(Ptr::in, X::in, Success::out, Result::out,
        _IO0::di, _IO::uo),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    double result;
    int code = gelfgren_pade_evaluate((void *)Ptr, X, &result);
    Success = (code == 0) ? MR_YES : MR_NO;
    Result = result;
").

%---------------------------------------------------------------------------%
% Helper predicates
%---------------------------------------------------------------------------%

:- pred is_null_ptr(T::in) is semidet.

:- pragma foreign_proc("C",
    is_null_ptr(Ptr::in),
    [will_not_call_mercury, promise_pure, thread_safe],
"
    SUCCESS_INDICATOR = (Ptr == (MR_Word)NULL);
").

%---------------------------------------------------------------------------%
:- end_module gelfgren.
%---------------------------------------------------------------------------%
