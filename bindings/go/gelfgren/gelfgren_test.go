package gelfgren

import (
	"math"
	"testing"
)

const epsilon = 1e-10

func TestBernsteinPolynomialCreate(t *testing.T) {
	poly, err := NewBernsteinPolynomial([]float64{1.0, 2.0, 6.0}, 0.0, 1.0)
	if err != nil {
		t.Fatalf("Failed to create polynomial: %v", err)
	}
	defer poly.Free()

	deg, err := poly.Degree()
	if err != nil {
		t.Fatalf("Failed to get degree: %v", err)
	}
	if deg != 2 {
		t.Errorf("Expected degree 2, got %d", deg)
	}
}

func TestBernsteinPolynomialEvaluate(t *testing.T) {
	poly, err := NewBernsteinPolynomial([]float64{1.0, 2.0, 6.0}, 0.0, 1.0)
	if err != nil {
		t.Fatalf("Failed to create polynomial: %v", err)
	}
	defer poly.Free()

	y, err := poly.Evaluate(0.5)
	if err != nil {
		t.Fatalf("Failed to evaluate: %v", err)
	}
	if math.Abs(y-3.0) > epsilon {
		t.Errorf("Expected 3.0, got %f", y)
	}
}

func TestBernsteinPolynomialEvaluateSlice(t *testing.T) {
	poly, err := NewBernsteinPolynomial([]float64{1.0, 2.0, 6.0}, 0.0, 1.0)
	if err != nil {
		t.Fatalf("Failed to create polynomial: %v", err)
	}
	defer poly.Free()

	xs := []float64{0.0, 0.5, 1.0}
	ys, err := poly.EvaluateSlice(xs)
	if err != nil {
		t.Fatalf("Failed to evaluate slice: %v", err)
	}

	expected := []float64{1.0, 3.0, 6.0}
	for i := range xs {
		if math.Abs(ys[i]-expected[i]) > epsilon {
			t.Errorf("At x=%f: expected %f, got %f", xs[i], expected[i], ys[i])
		}
	}
}

func TestBernsteinPolynomialDerivative(t *testing.T) {
	poly, err := NewBernsteinPolynomial([]float64{1.0, 2.0, 6.0}, 0.0, 1.0)
	if err != nil {
		t.Fatalf("Failed to create polynomial: %v", err)
	}
	defer poly.Free()

	dpoly, err := poly.Derivative()
	if err != nil {
		t.Fatalf("Failed to compute derivative: %v", err)
	}
	defer dpoly.Free()

	dy, err := dpoly.Evaluate(0.5)
	if err != nil {
		t.Fatalf("Failed to evaluate derivative: %v", err)
	}

	// Derivative of 2x^2 + 3x + 1 is 4x + 3, at x=0.5 is 5.0
	if math.Abs(dy-5.0) > epsilon {
		t.Errorf("Expected 5.0, got %f", dy)
	}
}

func TestBernsteinPolynomialIntegral(t *testing.T) {
	poly, err := NewBernsteinPolynomial([]float64{1.0, 2.0, 6.0}, 0.0, 1.0)
	if err != nil {
		t.Fatalf("Failed to create polynomial: %v", err)
	}
	defer poly.Free()

	ipoly, err := poly.Integral()
	if err != nil {
		t.Fatalf("Failed to compute integral: %v", err)
	}
	defer ipoly.Free()

	y1, err := ipoly.Evaluate(1.0)
	if err != nil {
		t.Fatalf("Failed to evaluate integral at 1: %v", err)
	}

	y0, err := ipoly.Evaluate(0.0)
	if err != nil {
		t.Fatalf("Failed to evaluate integral at 0: %v", err)
	}

	area := y1 - y0
	expected := 2.0/3.0 + 3.0/2.0 + 1.0 // Integral of 2x^2 + 3x + 1 from 0 to 1
	if math.Abs(area-expected) > epsilon {
		t.Errorf("Expected %f, got %f", expected, area)
	}
}

func TestBernsteinPolynomialArithmetic(t *testing.T) {
	p1, err := NewBernsteinPolynomial([]float64{1.0, 2.0, 3.0}, 0.0, 1.0)
	if err != nil {
		t.Fatalf("Failed to create p1: %v", err)
	}
	defer p1.Free()

	p2, err := NewBernsteinPolynomial([]float64{2.0, 1.0, 1.0}, 0.0, 1.0)
	if err != nil {
		t.Fatalf("Failed to create p2: %v", err)
	}
	defer p2.Free()

	y1, _ := p1.Evaluate(0.5)
	y2, _ := p2.Evaluate(0.5)

	// Test addition
	sum, err := p1.Add(p2)
	if err != nil {
		t.Fatalf("Failed to add: %v", err)
	}
	defer sum.Free()

	ySum, _ := sum.Evaluate(0.5)
	if math.Abs(ySum-(y1+y2)) > epsilon {
		t.Errorf("Addition failed: expected %f, got %f", y1+y2, ySum)
	}

	// Test scaling
	scaled, err := p1.Scale(2.0)
	if err != nil {
		t.Fatalf("Failed to scale: %v", err)
	}
	defer scaled.Free()

	yScaled, _ := scaled.Evaluate(0.5)
	if math.Abs(yScaled-2.0*y1) > epsilon {
		t.Errorf("Scaling failed: expected %f, got %f", 2.0*y1, yScaled)
	}
}

func TestRationalFunction(t *testing.T) {
	num, err := NewBernsteinPolynomial([]float64{0.0, 1.0}, 0.0, 1.0)
	if err != nil {
		t.Fatalf("Failed to create numerator: %v", err)
	}
	defer num.Free()

	den, err := NewBernsteinPolynomial([]float64{1.0, 2.0}, 0.0, 1.0)
	if err != nil {
		t.Fatalf("Failed to create denominator: %v", err)
	}
	defer den.Free()

	rat, err := NewRationalFunction(num, den)
	if err != nil {
		t.Fatalf("Failed to create rational function: %v", err)
	}
	defer rat.Free()

	x := 0.5
	y, err := rat.Evaluate(x)
	if err != nil {
		t.Fatalf("Failed to evaluate rational: %v", err)
	}

	expected := x / (1.0 + x)
	if math.Abs(y-expected) > epsilon {
		t.Errorf("Expected %f, got %f", expected, y)
	}
}

func TestPadeApproximant(t *testing.T) {
	coeffs := []float64{1.0, 1.0, 0.5, 1.0 / 6.0, 1.0 / 24.0}
	pade, err := NewPadeApproximant(coeffs, 2, 2, -1.0, 1.0)
	if err != nil {
		t.Fatalf("Failed to create Padé approximant: %v", err)
	}
	defer pade.Free()

	x := 0.5
	y, err := pade.Evaluate(x)
	if err != nil {
		t.Fatalf("Failed to evaluate Padé: %v", err)
	}

	exact := math.Exp(x)
	error := math.Abs(y - exact)

	if error > 0.001 {
		t.Errorf("Padé approximation error too large: %f", error)
	}
}

func TestErrorHandling(t *testing.T) {
	_, err := NewBernsteinPolynomial([]float64{}, 0.0, 1.0)
	if err == nil {
		t.Error("Expected error for empty coefficients")
	}
}

func BenchmarkEvaluate(b *testing.B) {
	poly, _ := NewBernsteinPolynomial([]float64{1.0, 2.0, 6.0}, 0.0, 1.0)
	defer poly.Free()

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		poly.Evaluate(0.5)
	}
}

func BenchmarkEvaluateSlice(b *testing.B) {
	poly, _ := NewBernsteinPolynomial([]float64{1.0, 2.0, 6.0}, 0.0, 1.0)
	defer poly.Free()

	xs := make([]float64, 1000)
	for i := range xs {
		xs[i] = float64(i) / 1000.0
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		poly.EvaluateSlice(xs)
	}
}
