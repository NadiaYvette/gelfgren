! Test suite for Gelfgren Fortran bindings
program test_bernstein
    use gelfgren_mod
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    integer :: n_tests, n_passed, n_failed
    real(real64), parameter :: TOL = 1.0d-9

    n_tests = 0
    n_passed = 0
    n_failed = 0

    print *, "Gelfgren Fortran Bindings Test Suite"
    print *, repeat("=", 60)

    ! Run tests
    call test_bernstein_creation()
    call test_bernstein_evaluation()
    call test_bernstein_derivative()
    call test_bernstein_integral()
    call test_bernstein_arithmetic()
    call test_rational_creation()
    call test_rational_evaluation()
    call test_pade_creation()
    call test_pade_evaluation()

    ! Summary
    print *
    print *, repeat("=", 60)
    print '(A,I0,A,I0,A,I0)', "Tests: ", n_tests, " passed, ", n_passed, &
                              " failed, ", n_failed
    if (n_failed == 0) then
        print *, "All tests passed!"
    else
        print *, "Some tests failed!"
        stop 1
    end if

contains

    subroutine test_bernstein_creation()
        type(bernstein_polynomial) :: poly
        type(gelfgren_error) :: error
        real(real64) :: coeffs(3)

        print *
        print *, "Test: Bernstein polynomial creation"

        coeffs = [1.0d0, 2.0d0, 3.0d0]
        call poly%create(coeffs, 0.0d0, 1.0d0, error)

        n_tests = n_tests + 1
        if (error%code == GELFGREN_SUCCESS) then
            print *, "  ✓ Creation succeeded"
            n_passed = n_passed + 1
        else
            print *, "  ✗ Creation failed:", trim(error%message)
            n_failed = n_failed + 1
        end if

        call poly%free()
    end subroutine test_bernstein_creation

    subroutine test_bernstein_evaluation()
        type(bernstein_polynomial) :: poly
        type(gelfgren_error) :: error
        real(real64) :: coeffs(2), y

        print *
        print *, "Test: Bernstein polynomial evaluation"

        ! Constant polynomial: f(x) = 2
        coeffs = [2.0d0, 2.0d0]
        call poly%create(coeffs, 0.0d0, 1.0d0, error)

        y = poly%evaluate(0.5d0, error)

        n_tests = n_tests + 1
        if (abs(y - 2.0d0) < TOL) then
            print *, "  ✓ Evaluation correct: f(0.5) = ", y
            n_passed = n_passed + 1
        else
            print *, "  ✗ Evaluation incorrect: expected 2.0, got", y
            n_failed = n_failed + 1
        end if

        ! Test endpoint interpolation
        y = poly%evaluate(0.0d0, error)
        n_tests = n_tests + 1
        if (abs(y - 2.0d0) < TOL) then
            print *, "  ✓ Left endpoint: f(0) = ", y
            n_passed = n_passed + 1
        else
            print *, "  ✗ Left endpoint incorrect"
            n_failed = n_failed + 1
        end if

        y = poly%evaluate(1.0d0, error)
        n_tests = n_tests + 1
        if (abs(y - 2.0d0) < TOL) then
            print *, "  ✓ Right endpoint: f(1) = ", y
            n_passed = n_passed + 1
        else
            print *, "  ✗ Right endpoint incorrect"
            n_failed = n_failed + 1
        end if

        call poly%free()
    end subroutine test_bernstein_evaluation

    subroutine test_bernstein_derivative()
        type(bernstein_polynomial) :: poly, dpoly
        type(gelfgren_error) :: error
        real(real64) :: coeffs(2), y

        print *
        print *, "Test: Bernstein polynomial derivative"

        ! Linear: f(x) = 2x + 1, f'(x) = 2
        coeffs = [1.0d0, 3.0d0]
        call poly%create(coeffs, 0.0d0, 1.0d0, error)

        dpoly = poly%derivative()
        y = dpoly%evaluate(0.5d0, error)

        n_tests = n_tests + 1
        if (abs(y - 2.0d0) < TOL) then
            print *, "  ✓ Derivative correct: f'(0.5) = ", y
            n_passed = n_passed + 1
        else
            print *, "  ✗ Derivative incorrect: expected 2.0, got", y
            n_failed = n_failed + 1
        end if

        call poly%free()
        call dpoly%free()
    end subroutine test_bernstein_derivative

    subroutine test_bernstein_integral()
        type(bernstein_polynomial) :: poly, ipoly
        type(gelfgren_error) :: error
        real(real64) :: coeffs(2), area

        print *
        print *, "Test: Bernstein polynomial integral"

        ! Constant: f(x) = 2, ∫f = 2x
        coeffs = [2.0d0, 2.0d0]
        call poly%create(coeffs, 0.0d0, 1.0d0, error)

        ipoly = poly%integral()
        area = ipoly%evaluate(1.0d0, error) - ipoly%evaluate(0.0d0, error)

        n_tests = n_tests + 1
        if (abs(area - 2.0d0) < TOL) then
            print *, "  ✓ Integral correct: area = ", area
            n_passed = n_passed + 1
        else
            print *, "  ✗ Integral incorrect: expected 2.0, got", area
            n_failed = n_failed + 1
        end if

        call poly%free()
        call ipoly%free()
    end subroutine test_bernstein_integral

    subroutine test_bernstein_arithmetic()
        type(bernstein_polynomial) :: p1, p2, result_poly
        type(gelfgren_error) :: error
        real(real64) :: coeffs(2), y

        print *
        print *, "Test: Bernstein polynomial arithmetic"

        coeffs = [1.0d0, 1.0d0]
        call p1%create(coeffs, 0.0d0, 1.0d0, error)

        coeffs = [2.0d0, 2.0d0]
        call p2%create(coeffs, 0.0d0, 1.0d0, error)

        ! Addition: 1 + 2 = 3
        result_poly = p1%add(p2)
        y = result_poly%evaluate(0.5d0, error)

        n_tests = n_tests + 1
        if (abs(y - 3.0d0) < TOL) then
            print *, "  ✓ Addition correct: p1 + p2 = ", y
            n_passed = n_passed + 1
        else
            print *, "  ✗ Addition incorrect: expected 3.0, got", y
            n_failed = n_failed + 1
        end if
        call result_poly%free()

        ! Subtraction: 2 - 1 = 1
        result_poly = p2%subtract(p1)
        y = result_poly%evaluate(0.5d0, error)

        n_tests = n_tests + 1
        if (abs(y - 1.0d0) < TOL) then
            print *, "  ✓ Subtraction correct: p2 - p1 = ", y
            n_passed = n_passed + 1
        else
            print *, "  ✗ Subtraction incorrect: expected 1.0, got", y
            n_failed = n_failed + 1
        end if
        call result_poly%free()

        ! Scaling: 2 * 2 = 4
        result_poly = p2%scale(2.0d0)
        y = result_poly%evaluate(0.5d0, error)

        n_tests = n_tests + 1
        if (abs(y - 4.0d0) < TOL) then
            print *, "  ✓ Scaling correct: 2 * p2 = ", y
            n_passed = n_passed + 1
        else
            print *, "  ✗ Scaling incorrect: expected 4.0, got", y
            n_failed = n_failed + 1
        end if
        call result_poly%free()

        call p1%free()
        call p2%free()
    end subroutine test_bernstein_arithmetic

    subroutine test_rational_creation()
        type(bernstein_polynomial) :: num, den
        type(rational_function) :: rat
        type(gelfgren_error) :: error
        real(real64) :: coeffs(2)

        print *
        print *, "Test: Rational function creation"

        coeffs = [0.0d0, 1.0d0]
        call num%create(coeffs, 0.0d0, 1.0d0, error)

        coeffs = [1.0d0, 1.0d0]
        call den%create(coeffs, 0.0d0, 1.0d0, error)

        call rat%create(num, den, error)

        n_tests = n_tests + 1
        if (error%code == GELFGREN_SUCCESS) then
            print *, "  ✓ Creation succeeded"
            n_passed = n_passed + 1
        else
            print *, "  ✗ Creation failed:", trim(error%message)
            n_failed = n_failed + 1
        end if

        call num%free()
        call den%free()
        call rat%free()
    end subroutine test_rational_creation

    subroutine test_rational_evaluation()
        type(bernstein_polynomial) :: num, den
        type(rational_function) :: rat
        type(gelfgren_error) :: error
        real(real64) :: coeffs(2), y, expected

        print *
        print *, "Test: Rational function evaluation"

        ! x / 1 = x
        coeffs = [0.0d0, 1.0d0]
        call num%create(coeffs, 0.0d0, 1.0d0, error)

        coeffs = [1.0d0, 1.0d0]
        call den%create(coeffs, 0.0d0, 1.0d0, error)

        call rat%create(num, den, error)

        y = rat%evaluate(0.5d0, error)
        expected = 0.5d0

        n_tests = n_tests + 1
        if (abs(y - expected) < TOL) then
            print *, "  ✓ Evaluation correct: R(0.5) = ", y
            n_passed = n_passed + 1
        else
            print *, "  ✗ Evaluation incorrect: expected", expected, ", got", y
            n_failed = n_failed + 1
        end if

        call num%free()
        call den%free()
        call rat%free()
    end subroutine test_rational_evaluation

    subroutine test_pade_creation()
        type(pade_approximant) :: pade
        type(gelfgren_error) :: error
        real(real64) :: coeffs(5)

        print *
        print *, "Test: Padé approximant creation"

        coeffs = [1.0d0, 1.0d0, 0.5d0, 1.0d0/6.0d0, 1.0d0/24.0d0]
        call pade%create_from_series(coeffs, 2, 2, -1.0d0, 1.0d0, error)

        n_tests = n_tests + 1
        if (error%code == GELFGREN_SUCCESS) then
            print *, "  ✓ Creation succeeded"
            n_passed = n_passed + 1
        else
            print *, "  ✗ Creation failed:", trim(error%message)
            n_failed = n_failed + 1
        end if

        call pade%free()
    end subroutine test_pade_creation

    subroutine test_pade_evaluation()
        type(pade_approximant) :: pade
        type(gelfgren_error) :: error
        real(real64) :: coeffs(5), y, exact, error_val

        print *
        print *, "Test: Padé approximant evaluation"

        ! Approximate exp(x) with [2/2]
        coeffs = [1.0d0, 1.0d0, 0.5d0, 1.0d0/6.0d0, 1.0d0/24.0d0]
        call pade%create_from_series(coeffs, 2, 2, -1.0d0, 1.0d0, error)

        y = pade%evaluate(0.5d0, error)
        exact = exp(0.5d0)
        error_val = abs(y - exact)

        n_tests = n_tests + 1
        if (error_val < 1.0d-5) then
            print '(A,E12.4)', "  ✓ Approximation accurate (error: ", error_val, ")"
            n_passed = n_passed + 1
        else
            print *, "  ✗ Approximation inaccurate (error:", error_val, ")"
            n_failed = n_failed + 1
        end if

        call pade%free()
    end subroutine test_pade_evaluation

end program test_bernstein
