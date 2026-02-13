! Gelfgren Fortran Bindings
! High-level Fortran interface to the Gelfgren C API
module gelfgren_mod
    use, intrinsic :: iso_c_binding
    implicit none

    private
    public :: bernstein_polynomial, rational_function, pade_approximant
    public :: gelfgren_error

    ! Error codes
    integer, parameter, public :: GELFGREN_SUCCESS = 0
    integer, parameter, public :: GELFGREN_NULL_POINTER = -1
    integer, parameter, public :: GELFGREN_INVALID_INTERVAL = -2
    integer, parameter, public :: GELFGREN_EMPTY_COEFFICIENTS = -3
    integer, parameter, public :: GELFGREN_SINGULAR_MATRIX = -4
    integer, parameter, public :: GELFGREN_DIVISION_BY_ZERO = -5
    integer, parameter, public :: GELFGREN_POLE = -6

    ! Opaque pointer types
    type, public :: bernstein_polynomial
        type(c_ptr) :: ptr = c_null_ptr
    contains
        procedure :: create => bernstein_create
        procedure :: free => bernstein_free
        procedure :: evaluate_scalar => bernstein_evaluate_scalar
        procedure :: evaluate_array => bernstein_evaluate_array
        generic :: evaluate => evaluate_scalar, evaluate_array
        procedure :: derivative => bernstein_derivative
        procedure :: integral => bernstein_integral
        procedure :: degree => bernstein_degree
        procedure :: add => bernstein_add
        procedure :: subtract => bernstein_subtract
        procedure :: multiply => bernstein_multiply
        procedure :: scale => bernstein_scale
    end type bernstein_polynomial

    type, public :: rational_function
        type(c_ptr) :: ptr = c_null_ptr
    contains
        procedure :: create => rational_create
        procedure :: free => rational_free
        procedure :: evaluate_scalar => rational_evaluate_scalar
        procedure :: evaluate_array => rational_evaluate_array
        generic :: evaluate => evaluate_scalar, evaluate_array
        procedure :: derivative => rational_derivative
    end type rational_function

    type, public :: pade_approximant
        type(c_ptr) :: ptr = c_null_ptr
    contains
        procedure :: create_from_series => pade_create_from_series
        procedure :: free => pade_free
        procedure :: evaluate_scalar => pade_evaluate_scalar
        procedure :: evaluate_array => pade_evaluate_array
        generic :: evaluate => evaluate_scalar, evaluate_array
    end type pade_approximant

    ! Error handling
    type, public :: gelfgren_error
        integer :: code
        character(len=256) :: message
    end type gelfgren_error

    ! C interface declarations
    interface
        ! Bernstein polynomial functions
        function gelfgren_bernstein_create(coeffs, degree, a, b) bind(C, name='gelfgren_bernstein_create')
            import :: c_ptr, c_double, c_size_t
            real(c_double), intent(in) :: coeffs(*)
            integer(c_size_t), value :: degree
            real(c_double), value :: a, b
            type(c_ptr) :: gelfgren_bernstein_create
        end function gelfgren_bernstein_create

        subroutine gelfgren_bernstein_free(ptr) bind(C, name='gelfgren_bernstein_free')
            import :: c_ptr
            type(c_ptr), value :: ptr
        end subroutine gelfgren_bernstein_free

        function gelfgren_bernstein_evaluate(ptr, x, result) bind(C, name='gelfgren_bernstein_evaluate')
            import :: c_ptr, c_double, c_int
            type(c_ptr), value :: ptr
            real(c_double), value :: x
            real(c_double) :: result
            integer(c_int) :: gelfgren_bernstein_evaluate
        end function gelfgren_bernstein_evaluate

        function gelfgren_bernstein_derivative(ptr) bind(C, name='gelfgren_bernstein_derivative')
            import :: c_ptr
            type(c_ptr), value :: ptr
            type(c_ptr) :: gelfgren_bernstein_derivative
        end function gelfgren_bernstein_derivative

        function gelfgren_bernstein_integral(ptr) bind(C, name='gelfgren_bernstein_integral')
            import :: c_ptr
            type(c_ptr), value :: ptr
            type(c_ptr) :: gelfgren_bernstein_integral
        end function gelfgren_bernstein_integral

        function gelfgren_bernstein_degree(ptr, degree) bind(C, name='gelfgren_bernstein_degree')
            import :: c_ptr, c_size_t, c_int
            type(c_ptr), value :: ptr
            integer(c_size_t) :: degree
            integer(c_int) :: gelfgren_bernstein_degree
        end function gelfgren_bernstein_degree

        function gelfgren_bernstein_add(ptr1, ptr2) bind(C, name='gelfgren_bernstein_add')
            import :: c_ptr
            type(c_ptr), value :: ptr1, ptr2
            type(c_ptr) :: gelfgren_bernstein_add
        end function gelfgren_bernstein_add

        function gelfgren_bernstein_subtract(ptr1, ptr2) bind(C, name='gelfgren_bernstein_subtract')
            import :: c_ptr
            type(c_ptr), value :: ptr1, ptr2
            type(c_ptr) :: gelfgren_bernstein_subtract
        end function gelfgren_bernstein_subtract

        function gelfgren_bernstein_multiply(ptr1, ptr2) bind(C, name='gelfgren_bernstein_multiply')
            import :: c_ptr
            type(c_ptr), value :: ptr1, ptr2
            type(c_ptr) :: gelfgren_bernstein_multiply
        end function gelfgren_bernstein_multiply

        function gelfgren_bernstein_scale(ptr, scalar) bind(C, name='gelfgren_bernstein_scale')
            import :: c_ptr, c_double
            type(c_ptr), value :: ptr
            real(c_double), value :: scalar
            type(c_ptr) :: gelfgren_bernstein_scale
        end function gelfgren_bernstein_scale

        ! Rational function functions
        function gelfgren_rational_create(num, den) bind(C, name='gelfgren_rational_create')
            import :: c_ptr
            type(c_ptr), value :: num, den
            type(c_ptr) :: gelfgren_rational_create
        end function gelfgren_rational_create

        subroutine gelfgren_rational_free(ptr) bind(C, name='gelfgren_rational_free')
            import :: c_ptr
            type(c_ptr), value :: ptr
        end subroutine gelfgren_rational_free

        function gelfgren_rational_evaluate(ptr, x, result) bind(C, name='gelfgren_rational_evaluate')
            import :: c_ptr, c_double, c_int
            type(c_ptr), value :: ptr
            real(c_double), value :: x
            real(c_double) :: result
            integer(c_int) :: gelfgren_rational_evaluate
        end function gelfgren_rational_evaluate

        function gelfgren_rational_derivative(ptr) bind(C, name='gelfgren_rational_derivative')
            import :: c_ptr
            type(c_ptr), value :: ptr
            type(c_ptr) :: gelfgren_rational_derivative
        end function gelfgren_rational_derivative

        ! Padé approximant functions
        function gelfgren_pade_from_power_series(coeffs, n_coeffs, n, m, a, b) &
                bind(C, name='gelfgren_pade_from_power_series')
            import :: c_ptr, c_double, c_size_t
            real(c_double), intent(in) :: coeffs(*)
            integer(c_size_t), value :: n_coeffs, n, m
            real(c_double), value :: a, b
            type(c_ptr) :: gelfgren_pade_from_power_series
        end function gelfgren_pade_from_power_series

        subroutine gelfgren_pade_free(ptr) bind(C, name='gelfgren_pade_free')
            import :: c_ptr
            type(c_ptr), value :: ptr
        end subroutine gelfgren_pade_free

        function gelfgren_pade_evaluate(ptr, x, result) bind(C, name='gelfgren_pade_evaluate')
            import :: c_ptr, c_double, c_int
            type(c_ptr), value :: ptr
            real(c_double), value :: x
            real(c_double) :: result
            integer(c_int) :: gelfgren_pade_evaluate
        end function gelfgren_pade_evaluate

        ! Error handling
        function gelfgren_last_error_message() bind(C, name='gelfgren_last_error_message')
            import :: c_ptr
            type(c_ptr) :: gelfgren_last_error_message
        end function gelfgren_last_error_message
    end interface

contains

    ! Bernstein polynomial methods
    subroutine bernstein_create(self, coeffs, a, b, error)
        class(bernstein_polynomial), intent(inout) :: self
        real(c_double), intent(in) :: coeffs(:)
        real(c_double), intent(in) :: a, b
        type(gelfgren_error), intent(out), optional :: error

        integer(c_size_t) :: degree

        degree = size(coeffs, kind=c_size_t) - 1
        self%ptr = gelfgren_bernstein_create(coeffs, degree, a, b)

        if (.not. c_associated(self%ptr)) then
            if (present(error)) then
                error%code = GELFGREN_NULL_POINTER
                error%message = get_last_error_message()
            end if
        else
            if (present(error)) then
                error%code = GELFGREN_SUCCESS
                error%message = ""
            end if
        end if
    end subroutine bernstein_create

    subroutine bernstein_free(self)
        class(bernstein_polynomial), intent(inout) :: self
        if (c_associated(self%ptr)) then
            call gelfgren_bernstein_free(self%ptr)
            self%ptr = c_null_ptr
        end if
    end subroutine bernstein_free

    function bernstein_evaluate_scalar(self, x, error) result(y)
        class(bernstein_polynomial), intent(in) :: self
        real(c_double), intent(in) :: x
        type(gelfgren_error), intent(out), optional :: error
        real(c_double) :: y

        integer(c_int) :: err_code

        err_code = gelfgren_bernstein_evaluate(self%ptr, x, y)

        if (present(error)) then
            error%code = err_code
            if (err_code /= GELFGREN_SUCCESS) then
                error%message = get_last_error_message()
            else
                error%message = ""
            end if
        end if
    end function bernstein_evaluate_scalar

    function bernstein_evaluate_array(self, x, error) result(y)
        class(bernstein_polynomial), intent(in) :: self
        real(c_double), intent(in) :: x(:)
        type(gelfgren_error), intent(out), optional :: error
        real(c_double), allocatable :: y(:)

        integer :: i, n
        integer(c_int) :: err_code

        n = size(x)
        allocate(y(n))

        do i = 1, n
            err_code = gelfgren_bernstein_evaluate(self%ptr, x(i), y(i))
            if (err_code /= GELFGREN_SUCCESS) then
                if (present(error)) then
                    error%code = err_code
                    error%message = get_last_error_message()
                end if
                return
            end if
        end do

        if (present(error)) then
            error%code = GELFGREN_SUCCESS
            error%message = ""
        end if
    end function bernstein_evaluate_array

    function bernstein_derivative(self) result(dpoly)
        class(bernstein_polynomial), intent(in) :: self
        type(bernstein_polynomial) :: dpoly

        dpoly%ptr = gelfgren_bernstein_derivative(self%ptr)
    end function bernstein_derivative

    function bernstein_integral(self) result(ipoly)
        class(bernstein_polynomial), intent(in) :: self
        type(bernstein_polynomial) :: ipoly

        ipoly%ptr = gelfgren_bernstein_integral(self%ptr)
    end function bernstein_integral

    function bernstein_degree(self, error) result(deg)
        class(bernstein_polynomial), intent(in) :: self
        type(gelfgren_error), intent(out), optional :: error
        integer :: deg

        integer(c_size_t) :: degree
        integer(c_int) :: err_code

        err_code = gelfgren_bernstein_degree(self%ptr, degree)
        deg = int(degree)

        if (present(error)) then
            error%code = err_code
            if (err_code /= GELFGREN_SUCCESS) then
                error%message = get_last_error_message()
            else
                error%message = ""
            end if
        end if
    end function bernstein_degree

    function bernstein_add(self, other) result(sum_poly)
        class(bernstein_polynomial), intent(in) :: self, other
        type(bernstein_polynomial) :: sum_poly

        sum_poly%ptr = gelfgren_bernstein_add(self%ptr, other%ptr)
    end function bernstein_add

    function bernstein_subtract(self, other) result(diff_poly)
        class(bernstein_polynomial), intent(in) :: self, other
        type(bernstein_polynomial) :: diff_poly

        diff_poly%ptr = gelfgren_bernstein_subtract(self%ptr, other%ptr)
    end function bernstein_subtract

    function bernstein_multiply(self, other) result(prod_poly)
        class(bernstein_polynomial), intent(in) :: self, other
        type(bernstein_polynomial) :: prod_poly

        prod_poly%ptr = gelfgren_bernstein_multiply(self%ptr, other%ptr)
    end function bernstein_multiply

    function bernstein_scale(self, scalar) result(scaled_poly)
        class(bernstein_polynomial), intent(in) :: self
        real(c_double), intent(in) :: scalar
        type(bernstein_polynomial) :: scaled_poly

        scaled_poly%ptr = gelfgren_bernstein_scale(self%ptr, scalar)
    end function bernstein_scale

    ! Rational function methods
    subroutine rational_create(self, numerator, denominator, error)
        class(rational_function), intent(inout) :: self
        type(bernstein_polynomial), intent(in) :: numerator, denominator
        type(gelfgren_error), intent(out), optional :: error

        self%ptr = gelfgren_rational_create(numerator%ptr, denominator%ptr)

        if (.not. c_associated(self%ptr)) then
            if (present(error)) then
                error%code = GELFGREN_NULL_POINTER
                error%message = get_last_error_message()
            end if
        else
            if (present(error)) then
                error%code = GELFGREN_SUCCESS
                error%message = ""
            end if
        end if
    end subroutine rational_create

    subroutine rational_free(self)
        class(rational_function), intent(inout) :: self
        if (c_associated(self%ptr)) then
            call gelfgren_rational_free(self%ptr)
            self%ptr = c_null_ptr
        end if
    end subroutine rational_free

    function rational_evaluate_scalar(self, x, error) result(y)
        class(rational_function), intent(in) :: self
        real(c_double), intent(in) :: x
        type(gelfgren_error), intent(out), optional :: error
        real(c_double) :: y

        integer(c_int) :: err_code

        err_code = gelfgren_rational_evaluate(self%ptr, x, y)

        if (present(error)) then
            error%code = err_code
            if (err_code /= GELFGREN_SUCCESS) then
                error%message = get_last_error_message()
            else
                error%message = ""
            end if
        end if
    end function rational_evaluate_scalar

    function rational_evaluate_array(self, x, error) result(y)
        class(rational_function), intent(in) :: self
        real(c_double), intent(in) :: x(:)
        type(gelfgren_error), intent(out), optional :: error
        real(c_double), allocatable :: y(:)

        integer :: i, n
        integer(c_int) :: err_code

        n = size(x)
        allocate(y(n))

        do i = 1, n
            err_code = gelfgren_rational_evaluate(self%ptr, x(i), y(i))
            if (err_code /= GELFGREN_SUCCESS) then
                if (present(error)) then
                    error%code = err_code
                    error%message = get_last_error_message()
                end if
                return
            end if
        end do

        if (present(error)) then
            error%code = GELFGREN_SUCCESS
            error%message = ""
        end if
    end function rational_evaluate_array

    function rational_derivative(self) result(drat)
        class(rational_function), intent(in) :: self
        type(rational_function) :: drat

        drat%ptr = gelfgren_rational_derivative(self%ptr)
    end function rational_derivative

    ! Padé approximant methods
    subroutine pade_create_from_series(self, coeffs, n, m, a, b, error)
        class(pade_approximant), intent(inout) :: self
        real(c_double), intent(in) :: coeffs(:)
        integer, intent(in) :: n, m
        real(c_double), intent(in) :: a, b
        type(gelfgren_error), intent(out), optional :: error

        integer(c_size_t) :: n_coeffs

        n_coeffs = size(coeffs, kind=c_size_t)
        self%ptr = gelfgren_pade_from_power_series(coeffs, n_coeffs, &
                                                     int(n, c_size_t), int(m, c_size_t), a, b)

        if (.not. c_associated(self%ptr)) then
            if (present(error)) then
                error%code = GELFGREN_NULL_POINTER
                error%message = get_last_error_message()
            end if
        else
            if (present(error)) then
                error%code = GELFGREN_SUCCESS
                error%message = ""
            end if
        end if
    end subroutine pade_create_from_series

    subroutine pade_free(self)
        class(pade_approximant), intent(inout) :: self
        if (c_associated(self%ptr)) then
            call gelfgren_pade_free(self%ptr)
            self%ptr = c_null_ptr
        end if
    end subroutine pade_free

    function pade_evaluate_scalar(self, x, error) result(y)
        class(pade_approximant), intent(in) :: self
        real(c_double), intent(in) :: x
        type(gelfgren_error), intent(out), optional :: error
        real(c_double) :: y

        integer(c_int) :: err_code

        err_code = gelfgren_pade_evaluate(self%ptr, x, y)

        if (present(error)) then
            error%code = err_code
            if (err_code /= GELFGREN_SUCCESS) then
                error%message = get_last_error_message()
            else
                error%message = ""
            end if
        end if
    end function pade_evaluate_scalar

    function pade_evaluate_array(self, x, error) result(y)
        class(pade_approximant), intent(in) :: self
        real(c_double), intent(in) :: x(:)
        type(gelfgren_error), intent(out), optional :: error
        real(c_double), allocatable :: y(:)

        integer :: i, n
        integer(c_int) :: err_code

        n = size(x)
        allocate(y(n))

        do i = 1, n
            err_code = gelfgren_pade_evaluate(self%ptr, x(i), y(i))
            if (err_code /= GELFGREN_SUCCESS) then
                if (present(error)) then
                    error%code = err_code
                    error%message = get_last_error_message()
                end if
                return
            end if
        end do

        if (present(error)) then
            error%code = GELFGREN_SUCCESS
            error%message = ""
        end if
    end function pade_evaluate_array

    ! Helper functions
    function get_last_error_message() result(msg)
        character(len=256) :: msg
        type(c_ptr) :: c_msg_ptr
        character(kind=c_char), dimension(:), pointer :: c_msg
        integer :: i, msg_len

        c_msg_ptr = gelfgren_last_error_message()

        if (c_associated(c_msg_ptr)) then
            call c_f_pointer(c_msg_ptr, c_msg, [256])
            msg = ""
            do i = 1, 256
                if (c_msg(i) == c_null_char) exit
                msg(i:i) = c_msg(i)
            end do
        else
            msg = "Unknown error"
        end if
    end function get_last_error_message

end module gelfgren_mod
