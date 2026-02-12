// Package gelfgren provides Go bindings for the Gelfgren piecewise rational interpolation library.
package gelfgren

/*
#cgo CFLAGS: -I../../../include
#cgo LDFLAGS: -L../../../target/release -lgelfgren
#include <stdlib.h>
#include <string.h>
#include "gelfgren.h"
*/
import "C"
import (
	"fmt"
	"runtime"
	"unsafe"
)

// Error represents a Gelfgren error
type Error struct {
	Message string
}

func (e *Error) Error() string {
	return e.Message
}

// getLastError retrieves the last error message from the C library
func getLastError() string {
	msg := C.gelfgren_last_error_message()
	if msg == nil {
		return "Unknown error"
	}
	return C.GoString(msg)
}

// checkError checks an error code and returns a Go error if non-zero
func checkError(code C.int, operation string) error {
	if code != 0 {
		return &Error{Message: fmt.Sprintf("%s: %s", operation, getLastError())}
	}
	return nil
}

// checkNull checks for null pointer and returns error if null
func checkNull(ptr unsafe.Pointer, operation string) error {
	if ptr == nil {
		return &Error{Message: fmt.Sprintf("%s: %s", operation, getLastError())}
	}
	return nil
}

//
// BernsteinPolynomial
//

// BernsteinPolynomial represents a Bernstein polynomial on interval [a, b]
type BernsteinPolynomial struct {
	ptr unsafe.Pointer
	a   float64
	b   float64
}

// NewBernsteinPolynomial creates a new Bernstein polynomial from coefficients
func NewBernsteinPolynomial(coeffs []float64, a, b float64) (*BernsteinPolynomial, error) {
	if len(coeffs) == 0 {
		return nil, &Error{Message: "empty coefficient slice"}
	}

	degree := len(coeffs) - 1
	coeffsPtr := (*C.double)(unsafe.Pointer(&coeffs[0]))

	ptr := C.gelfgren_bernstein_create(coeffsPtr, C.size_t(degree), C.double(a), C.double(b))
	if err := checkNull(ptr, "create_bernstein"); err != nil {
		return nil, err
	}

	poly := &BernsteinPolynomial{ptr: ptr, a: a, b: b}
	runtime.SetFinalizer(poly, (*BernsteinPolynomial).Free)
	return poly, nil
}

// Free releases the memory associated with the polynomial
func (p *BernsteinPolynomial) Free() {
	if p.ptr != nil {
		C.gelfgren_bernstein_free(p.ptr)
		p.ptr = nil
	}
}

// Evaluate evaluates the polynomial at point x
func (p *BernsteinPolynomial) Evaluate(x float64) (float64, error) {
	var result C.double
	code := C.gelfgren_bernstein_evaluate(p.ptr, C.double(x), &result)
	if err := checkError(code, "evaluate"); err != nil {
		return 0, err
	}
	return float64(result), nil
}

// EvaluateSlice evaluates the polynomial at multiple points (vectorized)
func (p *BernsteinPolynomial) EvaluateSlice(xs []float64) ([]float64, error) {
	results := make([]float64, len(xs))
	for i, x := range xs {
		y, err := p.Evaluate(x)
		if err != nil {
			return nil, err
		}
		results[i] = y
	}
	return results, nil
}

// Derivative computes the derivative of the polynomial
func (p *BernsteinPolynomial) Derivative() (*BernsteinPolynomial, error) {
	ptr := C.gelfgren_bernstein_derivative(p.ptr)
	if err := checkNull(ptr, "derivative"); err != nil {
		return nil, err
	}

	dpoly := &BernsteinPolynomial{ptr: ptr, a: p.a, b: p.b}
	runtime.SetFinalizer(dpoly, (*BernsteinPolynomial).Free)
	return dpoly, nil
}

// Integral computes the antiderivative of the polynomial
func (p *BernsteinPolynomial) Integral() (*BernsteinPolynomial, error) {
	ptr := C.gelfgren_bernstein_integral(p.ptr)
	if err := checkNull(ptr, "integral"); err != nil {
		return nil, err
	}

	ipoly := &BernsteinPolynomial{ptr: ptr, a: p.a, b: p.b}
	runtime.SetFinalizer(ipoly, (*BernsteinPolynomial).Free)
	return ipoly, nil
}

// Degree returns the degree of the polynomial
func (p *BernsteinPolynomial) Degree() (int, error) {
	var deg C.int
	code := C.gelfgren_bernstein_degree(p.ptr, &deg)
	if err := checkError(code, "degree"); err != nil {
		return 0, err
	}
	return int(deg), nil
}

// Add adds two polynomials
func (p *BernsteinPolynomial) Add(other *BernsteinPolynomial) (*BernsteinPolynomial, error) {
	ptr := C.gelfgren_bernstein_add(p.ptr, other.ptr)
	if err := checkNull(ptr, "add"); err != nil {
		return nil, err
	}

	result := &BernsteinPolynomial{ptr: ptr, a: p.a, b: p.b}
	runtime.SetFinalizer(result, (*BernsteinPolynomial).Free)
	return result, nil
}

// Subtract subtracts another polynomial from this one
func (p *BernsteinPolynomial) Subtract(other *BernsteinPolynomial) (*BernsteinPolynomial, error) {
	ptr := C.gelfgren_bernstein_subtract(p.ptr, other.ptr)
	if err := checkNull(ptr, "subtract"); err != nil {
		return nil, err
	}

	result := &BernsteinPolynomial{ptr: ptr, a: p.a, b: p.b}
	runtime.SetFinalizer(result, (*BernsteinPolynomial).Free)
	return result, nil
}

// Multiply multiplies two polynomials
func (p *BernsteinPolynomial) Multiply(other *BernsteinPolynomial) (*BernsteinPolynomial, error) {
	ptr := C.gelfgren_bernstein_multiply(p.ptr, other.ptr)
	if err := checkNull(ptr, "multiply"); err != nil {
		return nil, err
	}

	result := &BernsteinPolynomial{ptr: ptr, a: p.a, b: p.b}
	runtime.SetFinalizer(result, (*BernsteinPolynomial).Free)
	return result, nil
}

// Scale multiplies the polynomial by a scalar
func (p *BernsteinPolynomial) Scale(scalar float64) (*BernsteinPolynomial, error) {
	ptr := C.gelfgren_bernstein_scale(C.double(scalar), p.ptr)
	if err := checkNull(ptr, "scale"); err != nil {
		return nil, err
	}

	result := &BernsteinPolynomial{ptr: ptr, a: p.a, b: p.b}
	runtime.SetFinalizer(result, (*BernsteinPolynomial).Free)
	return result, nil
}

//
// RationalFunction
//

// RationalFunction represents a rational function P(x)/Q(x)
type RationalFunction struct {
	ptr unsafe.Pointer
}

// NewRationalFunction creates a new rational function from numerator and denominator
func NewRationalFunction(num, den *BernsteinPolynomial) (*RationalFunction, error) {
	ptr := C.gelfgren_rational_create(num.ptr, den.ptr)
	if err := checkNull(ptr, "create_rational"); err != nil {
		return nil, err
	}

	rat := &RationalFunction{ptr: ptr}
	runtime.SetFinalizer(rat, (*RationalFunction).Free)
	return rat, nil
}

// Free releases the memory associated with the rational function
func (r *RationalFunction) Free() {
	if r.ptr != nil {
		C.gelfgren_rational_free(r.ptr)
		r.ptr = nil
	}
}

// Evaluate evaluates the rational function at point x
func (r *RationalFunction) Evaluate(x float64) (float64, error) {
	var result C.double
	code := C.gelfgren_rational_evaluate(r.ptr, C.double(x), &result)
	if err := checkError(code, "evaluate_rational"); err != nil {
		return 0, err
	}
	return float64(result), nil
}

// EvaluateSlice evaluates the rational function at multiple points
func (r *RationalFunction) EvaluateSlice(xs []float64) ([]float64, error) {
	results := make([]float64, len(xs))
	for i, x := range xs {
		y, err := r.Evaluate(x)
		if err != nil {
			return nil, err
		}
		results[i] = y
	}
	return results, nil
}

// Derivative computes the derivative using the quotient rule
func (r *RationalFunction) Derivative() (*RationalFunction, error) {
	ptr := C.gelfgren_rational_derivative(r.ptr)
	if err := checkNull(ptr, "derivative_rational"); err != nil {
		return nil, err
	}

	drat := &RationalFunction{ptr: ptr}
	runtime.SetFinalizer(drat, (*RationalFunction).Free)
	return drat, nil
}

//
// PadeApproximant
//

// PadeApproximant represents a Padé approximant [n/m]
type PadeApproximant struct {
	ptr unsafe.Pointer
}

// NewPadeApproximant creates a Padé approximant from power series coefficients
func NewPadeApproximant(coeffs []float64, n, m int, a, b float64) (*PadeApproximant, error) {
	coeffsPtr := (*C.double)(unsafe.Pointer(&coeffs[0]))

	ptr := C.gelfgren_pade_from_series(
		coeffsPtr,
		C.size_t(len(coeffs)),
		C.int(n),
		C.int(m),
		C.double(a),
		C.double(b),
	)
	if err := checkNull(ptr, "create_pade"); err != nil {
		return nil, err
	}

	pade := &PadeApproximant{ptr: ptr}
	runtime.SetFinalizer(pade, (*PadeApproximant).Free)
	return pade, nil
}

// Free releases the memory associated with the Padé approximant
func (p *PadeApproximant) Free() {
	if p.ptr != nil {
		C.gelfgren_pade_free(p.ptr)
		p.ptr = nil
	}
}

// Evaluate evaluates the Padé approximant at point x
func (p *PadeApproximant) Evaluate(x float64) (float64, error) {
	var result C.double
	code := C.gelfgren_pade_evaluate(p.ptr, C.double(x), &result)
	if err := checkError(code, "evaluate_pade"); err != nil {
		return 0, err
	}
	return float64(result), nil
}

// EvaluateSlice evaluates the Padé approximant at multiple points
func (p *PadeApproximant) EvaluateSlice(xs []float64) ([]float64, error) {
	results := make([]float64, len(xs))
	for i, x := range xs {
		y, err := p.Evaluate(x)
		if err != nil {
			return nil, err
		}
		results[i] = y
	}
	return results, nil
}
