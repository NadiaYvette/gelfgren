package main

import (
	"fmt"
	"log"
	"math"

	"github.com/yourusername/gelfgren-go/gelfgren"
)

func main() {
	fmt.Println("Gelfgren Go Bindings Example")
	fmt.Println(string(make([]byte, 60, 60)[:60]))
	fmt.Println()

	// Example 1: Bernstein Polynomial
	fmt.Println("1. Bernstein Polynomial")
	fmt.Println("------------------------------------------------------------")

	poly, err := gelfgren.NewBernsteinPolynomial([]float64{1.0, 2.0, 6.0}, 0.0, 1.0)
	if err != nil {
		log.Fatalf("Error creating polynomial: %v", err)
	}
	defer poly.Free()

	deg, _ := poly.Degree()
	fmt.Printf("Created polynomial with degree: %d\n", deg)

	y, _ := poly.Evaluate(0.5)
	fmt.Printf("f(0.5) = %.6f\n", y)

	// Vectorized evaluation
	fmt.Println("\nVectorized evaluation:")
	xs := []float64{0.0, 0.25, 0.5, 0.75, 1.0}
	ys, _ := poly.EvaluateSlice(xs)
	for i, x := range xs {
		fmt.Printf("  f(%.2f) = %.6f\n", x, ys[i])
	}

	// Derivative
	fmt.Println()
	dpoly, _ := poly.Derivative()
	defer dpoly.Free()
	dy, _ := dpoly.Evaluate(0.5)
	fmt.Printf("Derivative at x=0.5: %.6f\n", dy)

	// Integral
	ipoly, _ := poly.Integral()
	defer ipoly.Free()
	y1, _ := ipoly.Evaluate(1.0)
	y0, _ := ipoly.Evaluate(0.0)
	area := y1 - y0
	fmt.Printf("Integral from 0 to 1: %.6f\n", area)

	fmt.Println()

	// Example 2: Polynomial Arithmetic
	fmt.Println("2. Polynomial Arithmetic")
	fmt.Println("------------------------------------------------------------")

	p1, _ := gelfgren.NewBernsteinPolynomial([]float64{1.0, 2.0, 3.0}, 0.0, 1.0)
	defer p1.Free()

	p2, _ := gelfgren.NewBernsteinPolynomial([]float64{2.0, 1.0, 1.0}, 0.0, 1.0)
	defer p2.Free()

	y1, _ := p1.Evaluate(0.5)
	y2, _ := p2.Evaluate(0.5)
	fmt.Printf("p1(0.5) = %.6f\n", y1)
	fmt.Printf("p2(0.5) = %.6f\n", y2)

	sum, _ := p1.Add(p2)
	defer sum.Free()
	ySum, _ := sum.Evaluate(0.5)
	fmt.Printf("p1 + p2 at x=0.5: %.6f\n", ySum)

	scaled, _ := p1.Scale(2.0)
	defer scaled.Free()
	yScaled, _ := scaled.Evaluate(0.5)
	fmt.Printf("2 * p1 at x=0.5: %.6f\n", yScaled)

	fmt.Println()

	// Example 3: Rational Function
	fmt.Println("3. Rational Function")
	fmt.Println("------------------------------------------------------------")

	num, _ := gelfgren.NewBernsteinPolynomial([]float64{0.0, 1.0}, 0.0, 1.0)
	defer num.Free()

	den, _ := gelfgren.NewBernsteinPolynomial([]float64{1.0, 2.0}, 0.0, 1.0)
	defer den.Free()

	rat, _ := gelfgren.NewRationalFunction(num, den)
	defer rat.Free()

	x := 0.5
	yRat, _ := rat.Evaluate(x)
	expected := x / (1.0 + x)
	fmt.Printf("R(%.1f) = %.6f\n", x, yRat)
	fmt.Printf("Expected: %.6f\n", expected)

	fmt.Println("\nVectorized evaluation:")
	xsRat := []float64{0.0, 0.25, 0.5, 0.75, 1.0}
	ysRat, _ := rat.EvaluateSlice(xsRat)
	for i, x := range xsRat {
		expected := x / (1.0 + x)
		fmt.Printf("  R(%.2f) = %.6f (expected: %.6f)\n", x, ysRat[i], expected)
	}

	fmt.Println()

	// Example 4: Padé Approximant
	fmt.Println("4. Padé Approximant")
	fmt.Println("------------------------------------------------------------")

	coeffs := []float64{1.0, 1.0, 0.5, 1.0 / 6.0, 1.0 / 24.0}
	pade, _ := gelfgren.NewPadeApproximant(coeffs, 2, 2, -1.0, 1.0)
	defer pade.Free()

	fmt.Println("Approximation of exp(x):")
	xsPade := []float64{-0.5, -0.25, 0.0, 0.25, 0.5}
	ysPade, _ := pade.EvaluateSlice(xsPade)

	for i, x := range xsPade {
		exact := math.Exp(x)
		error := math.Abs(ysPade[i] - exact)
		fmt.Printf("  x=%.2f: Padé=%.6f, exp(x)=%.6f, error=%.8f\n",
			x, ysPade[i], exact, error)
	}

	fmt.Println()
	fmt.Println("All examples completed successfully!")
}
