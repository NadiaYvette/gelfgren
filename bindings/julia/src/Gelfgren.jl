"""
    Gelfgren

Julia bindings for the Gelfgren piecewise rational interpolation library.

Provides high-performance numerical methods for constructing rational approximants
based on Jan Gelfgren's 1975 research on piecewise rational interpolation.
"""
module Gelfgren

export BernsteinPolynomial, RationalFunction, PadeApproximant
export evaluate, derivative, integral, degree
export GelfgrenError

# Load library
const libgelfgren = let
    lib_path = get(ENV, "GELFGREN_LIB", "")
    if !isempty(lib_path)
        lib_path
    else
        # Try common locations
        for dir in ["../../../target/release", "/usr/local/lib", "/usr/lib"]
            path = joinpath(@__DIR__, dir, "libgelfgren.so")
            if isfile(path)
                return path
            end
            # Also try .dylib on macOS
            path = joinpath(@__DIR__, dir, "libgelfgren.dylib")
            if isfile(path)
                return path
            end
        end
        # Fallback to system library
        "libgelfgren"
    end
end

# Exception type
struct GelfgrenError <: Exception
    msg::String
end

# Helper function to get last error message
function get_last_error()
    ptr = ccall((:gelfgren_last_error_message, libgelfgren), Ptr{UInt8}, ())
    if ptr == C_NULL
        return "Unknown error"
    end
    unsafe_string(ptr)
end

# Helper to check error codes
function check_error(code::Int32, operation::String)
    if code != 0
        error_msg = get_last_error()
        throw(GelfgrenError("$operation failed: $error_msg"))
    end
end

# Helper to check null pointers
function check_null(ptr::Ptr{T}, operation::String) where T
    if ptr == C_NULL
        error_msg = get_last_error()
        throw(GelfgrenError("$operation failed: $error_msg"))
    end
    ptr
end

#
# Bernstein Polynomial
#

"""
    BernsteinPolynomial

Represents a Bernstein polynomial on interval [a, b].

# Fields
- `ptr::Ptr{Cvoid}`: Pointer to the underlying C structure
- `a::Float64`: Left endpoint of interval
- `b::Float64`: Right endpoint of interval

# Example
```julia
poly = BernsteinPolynomial([1.0, 2.0, 6.0], 0.0, 1.0)
y = evaluate(poly, 0.5)
```
"""
mutable struct BernsteinPolynomial
    ptr::Ptr{Cvoid}
    a::Float64
    b::Float64

    function BernsteinPolynomial(coeffs::AbstractVector{Float64}, a::Float64, b::Float64)
        if isempty(coeffs)
            throw(GelfgrenError("Empty coefficient vector"))
        end

        degree = length(coeffs) - 1
        ptr = ccall(
            (:gelfgren_bernstein_create, libgelfgren),
            Ptr{Cvoid},
            (Ptr{Float64}, Csize_t, Float64, Float64),
            coeffs, degree, a, b
        )

        ptr = check_null(ptr, "create_bernstein")

        poly = new(ptr, a, b)
        finalizer(free!, poly)
        return poly
    end
end

function free!(poly::BernsteinPolynomial)
    if poly.ptr != C_NULL
        ccall((:gelfgren_bernstein_free, libgelfgren), Cvoid, (Ptr{Cvoid},), poly.ptr)
        poly.ptr = C_NULL
    end
end

"""
    evaluate(poly::BernsteinPolynomial, x::Real)

Evaluate the Bernstein polynomial at point x.
"""
function evaluate(poly::BernsteinPolynomial, x::Real)
    result = Ref{Float64}(0.0)
    code = ccall(
        (:gelfgren_bernstein_evaluate, libgelfgren),
        Int32,
        (Ptr{Cvoid}, Float64, Ptr{Float64}),
        poly.ptr, Float64(x), result
    )
    check_error(code, "evaluate_bernstein")
    return result[]
end

"""
    evaluate(poly::BernsteinPolynomial, xs::AbstractVector)

Evaluate the Bernstein polynomial at multiple points (vectorized).
"""
function evaluate(poly::BernsteinPolynomial, xs::AbstractVector{<:Real})
    return [evaluate(poly, x) for x in xs]
end

"""
    derivative(poly::BernsteinPolynomial)

Compute the derivative of the Bernstein polynomial.
Returns a new BernsteinPolynomial.
"""
function derivative(poly::BernsteinPolynomial)
    ptr = ccall(
        (:gelfgren_bernstein_derivative, libgelfgren),
        Ptr{Cvoid},
        (Ptr{Cvoid},),
        poly.ptr
    )
    ptr = check_null(ptr, "derivative")

    dpoly = BernsteinPolynomial.__new__(BernsteinPolynomial)
    dpoly.ptr = ptr
    dpoly.a = poly.a
    dpoly.b = poly.b
    finalizer(free!, dpoly)
    return dpoly
end

"""
    integral(poly::BernsteinPolynomial)

Compute the antiderivative of the Bernstein polynomial.
Returns a new BernsteinPolynomial.
"""
function integral(poly::BernsteinPolynomial)
    ptr = ccall(
        (:gelfgren_bernstein_integral, libgelfgren),
        Ptr{Cvoid},
        (Ptr{Cvoid},),
        poly.ptr
    )
    ptr = check_null(ptr, "integral")

    ipoly = BernsteinPolynomial.__new__(BernsteinPolynomial)
    ipoly.ptr = ptr
    ipoly.a = poly.a
    ipoly.b = poly.b
    finalizer(free!, ipoly)
    return ipoly
end

"""
    degree(poly::BernsteinPolynomial)

Get the degree of the Bernstein polynomial.
"""
function degree(poly::BernsteinPolynomial)
    deg = Ref{Int32}(0)
    code = ccall(
        (:gelfgren_bernstein_degree, libgelfgren),
        Int32,
        (Ptr{Cvoid}, Ptr{Int32}),
        poly.ptr, deg
    )
    check_error(code, "degree")
    return Int(deg[])
end

# Arithmetic operations
import Base: +, -, *, /

"""
    +(p1::BernsteinPolynomial, p2::BernsteinPolynomial)

Add two Bernstein polynomials.
"""
function Base.:+(p1::BernsteinPolynomial, p2::BernsteinPolynomial)
    ptr = ccall(
        (:gelfgren_bernstein_add, libgelfgren),
        Ptr{Cvoid},
        (Ptr{Cvoid}, Ptr{Cvoid}),
        p1.ptr, p2.ptr
    )
    ptr = check_null(ptr, "add")

    result = BernsteinPolynomial.__new__(BernsteinPolynomial)
    result.ptr = ptr
    result.a = p1.a
    result.b = p1.b
    finalizer(free!, result)
    return result
end

"""
    -(p1::BernsteinPolynomial, p2::BernsteinPolynomial)

Subtract two Bernstein polynomials.
"""
function Base.:-(p1::BernsteinPolynomial, p2::BernsteinPolynomial)
    ptr = ccall(
        (:gelfgren_bernstein_subtract, libgelfgren),
        Ptr{Cvoid},
        (Ptr{Cvoid}, Ptr{Cvoid}),
        p1.ptr, p2.ptr
    )
    ptr = check_null(ptr, "subtract")

    result = BernsteinPolynomial.__new__(BernsteinPolynomial)
    result.ptr = ptr
    result.a = p1.a
    result.b = p1.b
    finalizer(free!, result)
    return result
end

"""
    *(p1::BernsteinPolynomial, p2::BernsteinPolynomial)

Multiply two Bernstein polynomials.
"""
function Base.:*(p1::BernsteinPolynomial, p2::BernsteinPolynomial)
    ptr = ccall(
        (:gelfgren_bernstein_multiply, libgelfgren),
        Ptr{Cvoid},
        (Ptr{Cvoid}, Ptr{Cvoid}),
        p1.ptr, p2.ptr
    )
    ptr = check_null(ptr, "multiply")

    result = BernsteinPolynomial.__new__(BernsteinPolynomial)
    result.ptr = ptr
    result.a = p1.a
    result.b = p1.b
    finalizer(free!, result)
    return result
end

"""
    *(scalar::Real, poly::BernsteinPolynomial)

Scale a Bernstein polynomial by a scalar.
"""
function Base.:*(scalar::Real, poly::BernsteinPolynomial)
    ptr = ccall(
        (:gelfgren_bernstein_scale, libgelfgren),
        Ptr{Cvoid},
        (Float64, Ptr{Cvoid}),
        Float64(scalar), poly.ptr
    )
    ptr = check_null(ptr, "scale")

    result = BernsteinPolynomial.__new__(BernsteinPolynomial)
    result.ptr = ptr
    result.a = poly.a
    result.b = poly.b
    finalizer(free!, result)
    return result
end

Base.:*(poly::BernsteinPolynomial, scalar::Real) = scalar * poly

#
# Rational Function
#

"""
    RationalFunction

Represents a rational function P(x)/Q(x) where P and Q are Bernstein polynomials.
"""
mutable struct RationalFunction
    ptr::Ptr{Cvoid}

    function RationalFunction(num::BernsteinPolynomial, den::BernsteinPolynomial)
        ptr = ccall(
            (:gelfgren_rational_create, libgelfgren),
            Ptr{Cvoid},
            (Ptr{Cvoid}, Ptr{Cvoid}),
            num.ptr, den.ptr
        )
        ptr = check_null(ptr, "create_rational")

        rat = new(ptr)
        finalizer(free!, rat)
        return rat
    end
end

function free!(rat::RationalFunction)
    if rat.ptr != C_NULL
        ccall((:gelfgren_rational_free, libgelfgren), Cvoid, (Ptr{Cvoid},), rat.ptr)
        rat.ptr = C_NULL
    end
end

"""
    evaluate(rat::RationalFunction, x::Real)

Evaluate the rational function at point x.
"""
function evaluate(rat::RationalFunction, x::Real)
    result = Ref{Float64}(0.0)
    code = ccall(
        (:gelfgren_rational_evaluate, libgelfgren),
        Int32,
        (Ptr{Cvoid}, Float64, Ptr{Float64}),
        rat.ptr, Float64(x), result
    )
    check_error(code, "evaluate_rational")
    return result[]
end

"""
    evaluate(rat::RationalFunction, xs::AbstractVector)

Evaluate the rational function at multiple points (vectorized).
"""
function evaluate(rat::RationalFunction, xs::AbstractVector{<:Real})
    return [evaluate(rat, x) for x in xs]
end

"""
    derivative(rat::RationalFunction)

Compute the derivative of the rational function using the quotient rule.
"""
function derivative(rat::RationalFunction)
    ptr = ccall(
        (:gelfgren_rational_derivative, libgelfgren),
        Ptr{Cvoid},
        (Ptr{Cvoid},),
        rat.ptr
    )
    ptr = check_null(ptr, "derivative_rational")

    drat = RationalFunction.__new__(RationalFunction)
    drat.ptr = ptr
    finalizer(free!, drat)
    return drat
end

#
# Padé Approximant
#

"""
    PadeApproximant

Represents a Padé approximant [n/m] constructed from power series coefficients.
"""
mutable struct PadeApproximant
    ptr::Ptr{Cvoid}

    function PadeApproximant(coeffs::AbstractVector{Float64}, n::Int, m::Int, a::Float64, b::Float64)
        ptr = ccall(
            (:gelfgren_pade_from_series, libgelfgren),
            Ptr{Cvoid},
            (Ptr{Float64}, Csize_t, Int32, Int32, Float64, Float64),
            coeffs, length(coeffs), Int32(n), Int32(m), a, b
        )
        ptr = check_null(ptr, "create_pade")

        pade = new(ptr)
        finalizer(free!, pade)
        return pade
    end
end

function free!(pade::PadeApproximant)
    if pade.ptr != C_NULL
        ccall((:gelfgren_pade_free, libgelfgren), Cvoid, (Ptr{Cvoid},), pade.ptr)
        pade.ptr = C_NULL
    end
end

"""
    evaluate(pade::PadeApproximant, x::Real)

Evaluate the Padé approximant at point x.
"""
function evaluate(pade::PadeApproximant, x::Real)
    result = Ref{Float64}(0.0)
    code = ccall(
        (:gelfgren_pade_evaluate, libgelfgren),
        Int32,
        (Ptr{Cvoid}, Float64, Ptr{Float64}),
        pade.ptr, Float64(x), result
    )
    check_error(code, "evaluate_pade")
    return result[]
end

"""
    evaluate(pade::PadeApproximant, xs::AbstractVector)

Evaluate the Padé approximant at multiple points (vectorized).
"""
function evaluate(pade::PadeApproximant, xs::AbstractVector{<:Real})
    return [evaluate(pade, x) for x in xs]
end

# Custom constructors for internal use
Base.@pure BernsteinPolynomial.__new__(::Type{BernsteinPolynomial}) = ccall(:jl_new_struct_uninit, Any, (Any,), BernsteinPolynomial)::BernsteinPolynomial
Base.@pure RationalFunction.__new__(::Type{RationalFunction}) = ccall(:jl_new_struct_uninit, Any, (Any,), RationalFunction)::RationalFunction
Base.@pure PadeApproximant.__new__(::Type{PadeApproximant}) = ccall(:jl_new_struct_uninit, Any, (Any,), PadeApproximant)::PadeApproximant

end # module Gelfgren
