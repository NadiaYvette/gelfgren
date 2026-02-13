! Gelfgren Fortran Bindings Example
! Demonstrates basic usage of the Gelfgren library from Fortran
program basic_usage
    use gelfgren_mod
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    ! Variables
    type(bernstein_polynomial) :: poly, dpoly, ipoly, poly2, sum_poly
    type(rational_function) :: rat, drat
    type(pade_approximant) :: pade
    type(gelfgren_error) :: error
    real(real64), allocatable :: coeffs(:), x_array(:), y_array(:)
    real(real64) :: x, y, area
    integer :: i, deg

    print *, "Gelfgren Fortran Bindings Example"
    print *, repeat("=", 60)

    ! Example 1: Bernstein Polynomial
    print *
    print *, "1. Bernstein Polynomial"
    print *, repeat("-", 60)

    ! Create a quadratic polynomial: 2x^2 + 3x + 1 on [0, 1]
    ! Bernstein coefficients: [1, 2, 6]
    allocate(coeffs(3))
    coeffs = [1.0d0, 2.0d0, 6.0d0]

    call poly%create(coeffs, 0.0d0, 1.0d0, error)
    if (error%code /= GELFGREN_SUCCESS) then
        print *, "Error creating polynomial:", trim(error%message)
        stop 1
    end if
    print *, "Created polynomial with degree:", poly%degree(error)

    ! Evaluate at a point
    x = 0.5d0
    y = poly%evaluate(x, error)
    if (error%code == GELFGREN_SUCCESS) then
        print *, "f(0.5) =", y
    else
        print *, "Error evaluating:", trim(error%message)
    end if

    ! Vectorized evaluation
    allocate(x_array(5))
    x_array = [0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0]
    y_array = poly%evaluate(x_array, error)

    print *
    print *, "Vectorized evaluation:"
    do i = 1, size(x_array)
        print '(A,F6.3,A,F10.6)', "  f(", x_array(i), ") = ", y_array(i)
    end do

    ! Derivative
    dpoly = poly%derivative()
    y = dpoly%evaluate(0.5d0, error)
    print *
    print *, "Derivative at x=0.5:", y

    ! Integral
    ipoly = poly%integral()
    area = ipoly%evaluate(1.0d0, error) - ipoly%evaluate(0.0d0, error)
    print *, "Integral from 0 to 1:", area

    ! Clean up
    call dpoly%free()
    call ipoly%free()
    deallocate(x_array, y_array, coeffs)

    ! Example 2: Polynomial Arithmetic
    print *
    print *, "2. Polynomial Arithmetic"
    print *, repeat("-", 60)

    allocate(coeffs(3))
    coeffs = [1.0d0, 2.0d0, 3.0d0]
    call poly2%create(coeffs, 0.0d0, 1.0d0, error)

    print *, "p1(0.5) =", poly%evaluate(0.5d0, error)
    print *, "p2(0.5) =", poly2%evaluate(0.5d0, error)

    ! Addition
    sum_poly = poly%add(poly2)
    print *, "p1 + p2 at x=0.5:", sum_poly%evaluate(0.5d0, error)
    call sum_poly%free()

    ! Scaling
    sum_poly = poly%scale(2.0d0)
    print *, "2 * p1 at x=0.5:", sum_poly%evaluate(0.5d0, error)
    call sum_poly%free()

    call poly2%free()
    deallocate(coeffs)

    ! Example 3: Rational Function
    print *
    print *, "3. Rational Function"
    print *, repeat("-", 60)

    ! Create rational function: x / (1 + x) on [0, 1]
    type(bernstein_polynomial) :: num, den

    allocate(coeffs(2))
    coeffs = [0.0d0, 1.0d0]
    call num%create(coeffs, 0.0d0, 1.0d0, error)

    coeffs = [1.0d0, 2.0d0]
    call den%create(coeffs, 0.0d0, 1.0d0, error)

    call rat%create(num, den, error)
    if (error%code /= GELFGREN_SUCCESS) then
        print *, "Error creating rational function:", trim(error%message)
        stop 1
    end if

    ! Evaluate
    x = 0.5d0
    y = rat%evaluate(x, error)
    print *, "R(0.5) =", y
    print *, "Expected:", x / (1.0d0 + x)

    ! Vectorized evaluation
    deallocate(coeffs)
    allocate(x_array(5))
    x_array = [0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0]
    y_array = rat%evaluate(x_array, error)

    print *
    print *, "Vectorized evaluation:"
    do i = 1, size(x_array)
        print '(A,F6.3,A,F10.6,A,F10.6)', "  R(", x_array(i), ") = ", &
              y_array(i), " (expected: ", x_array(i)/(1.0d0+x_array(i)), ")"
    end do

    ! Derivative
    drat = rat%derivative()
    y = drat%evaluate(0.5d0, error)
    print *
    print *, "R'(0.5) =", y

    ! Clean up
    call num%free()
    call den%free()
    call rat%free()
    call drat%free()
    deallocate(x_array, y_array)

    ! Example 4: Padé Approximant
    print *
    print *, "4. Padé Approximant"
    print *, repeat("-", 60)

    ! Approximate exp(x) near x=0 with [2/2] Padé approximant
    ! Power series: 1 + x + x^2/2 + x^3/6 + x^4/24
    allocate(coeffs(5))
    coeffs = [1.0d0, 1.0d0, 0.5d0, 1.0d0/6.0d0, 1.0d0/24.0d0]

    call pade%create_from_series(coeffs, 2, 2, -1.0d0, 1.0d0, error)
    if (error%code /= GELFGREN_SUCCESS) then
        print *, "Error creating Padé approximant:", trim(error%message)
        stop 1
    end if

    print *, "Approximation of exp(x):"
    allocate(x_array(5))
    x_array = [-0.5d0, -0.25d0, 0.0d0, 0.25d0, 0.5d0]

    do i = 1, size(x_array)
        x = x_array(i)
        y = pade%evaluate(x, error)
        if (error%code == GELFGREN_SUCCESS) then
            print '(A,F6.2,A,F10.6,A,F10.6,A,E12.4)', &
                "  x=", x, ": Padé=", y, ", exp(x)=", exp(x), &
                ", error=", abs(y - exp(x))
        else
            print *, "  Error at x=", x, ":", trim(error%message)
        end if
    end do

    ! Clean up
    call pade%free()
    call poly%free()
    deallocate(coeffs, x_array)

    print *
    print *, "All examples completed successfully!"

end program basic_usage
