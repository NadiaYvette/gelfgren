-- | Example usage of Gelfgren Haskell bindings
module Main where

import Control.Exception (catch, SomeException)
import Text.Printf (printf)
import Gelfgren

main :: IO ()
main = do
    putStrLn "Gelfgren Haskell Bindings Example"
    putStrLn (replicate 60 '=')
    putStrLn ""

    example1
    example2
    example3
    example4
    example5

    putStrLn ""
    putStrLn "All examples completed successfully!"

-- Example 1: Bernstein Polynomial
example1 :: IO ()
example1 = do
    putStrLn "1. Bernstein Polynomial"
    putStrLn (replicate 60 '-')

    -- Create a quadratic polynomial: 2x^2 + 3x + 1 on [0, 1]
    -- Bernstein coefficients: [1, 2, 6]
    poly <- bernsteinPolynomial [1.0, 2.0, 6.0] 0.0 1.0

    -- Get degree
    d <- degree poly
    printf "Created polynomial with degree: %d\n" d

    -- Evaluate at a point
    y <- evaluate poly 0.5
    printf "f(0.5) = %.6f\n" y

    -- Vectorized evaluation
    let xs = [0.0, 0.25, 0.5, 0.75, 1.0]
    ys <- evaluateMany poly xs
    putStrLn "\nVectorized evaluation:"
    mapM_ (\(x, y') -> printf "  f(%.2f) = %.6f\n" x y') (zip xs ys)

    -- Derivative
    dpoly <- derivative poly
    dy <- evaluate dpoly 0.5
    printf "\nDerivative at x=0.5: %.6f\n" dy

    -- Integral
    ipoly <- integral poly
    y0 <- evaluate ipoly 0.0
    y1 <- evaluate ipoly 1.0
    let area = y1 - y0
    printf "Integral from 0 to 1: %.6f\n" area

    putStrLn ""

-- Example 2: Polynomial Arithmetic
example2 :: IO ()
example2 = do
    putStrLn "2. Polynomial Arithmetic"
    putStrLn (replicate 60 '-')

    poly1 <- bernsteinPolynomial [1.0, 2.0, 3.0] 0.0 1.0
    poly2 <- bernsteinPolynomial [2.0, 1.0, 1.0] 0.0 1.0

    y1 <- evaluate poly1 0.5
    y2 <- evaluate poly2 0.5
    printf "p1(0.5) = %.6f\n" y1
    printf "p2(0.5) = %.6f\n" y2

    -- Addition
    sumPoly <- add poly1 poly2
    ySum <- evaluate sumPoly 0.5
    printf "p1 + p2 at x=0.5: %.6f\n" ySum

    -- Subtraction
    diffPoly <- subtract poly1 poly2
    yDiff <- evaluate diffPoly 0.5
    printf "p1 - p2 at x=0.5: %.6f\n" yDiff

    -- Multiplication
    prodPoly <- multiply poly1 poly2
    yProd <- evaluate prodPoly 0.5
    printf "p1 * p2 at x=0.5: %.6f\n" yProd

    -- Scaling
    scaledPoly <- scale 2.0 poly1
    yScaled <- evaluate scaledPoly 0.5
    printf "2 * p1 at x=0.5: %.6f\n" yScaled

    putStrLn ""

-- Example 3: Rational Function
example3 :: IO ()
example3 = do
    putStrLn "3. Rational Function"
    putStrLn (replicate 60 '-')

    -- Create rational function: x / (1 + x) on [0, 1]
    num <- bernsteinPolynomial [0.0, 1.0] 0.0 1.0
    den <- bernsteinPolynomial [1.0, 2.0] 0.0 1.0
    rat <- rationalFunction num den

    -- Evaluate
    let x = 0.5
    y <- evaluate rat x
    let expected = x / (1.0 + x)
    printf "R(%.1f) = %.6f\n" x y
    printf "Expected: %.6f\n" expected

    -- Vectorized evaluation
    let xs = [0.0, 0.25, 0.5, 0.75, 1.0]
    ys <- evaluateMany rat xs
    let expecteds = map (\xi -> xi / (1.0 + xi)) xs
    putStrLn "\nVectorized evaluation:"
    mapM_ (\(xi, yi, exp') ->
        printf "  R(%.2f) = %.6f (expected: %.6f)\n" xi yi exp')
        (zip3 xs ys expecteds)

    putStrLn ""

-- Example 4: Padé Approximant
example4 :: IO ()
example4 = do
    putStrLn "4. Padé Approximant"
    putStrLn (replicate 60 '-')

    -- Approximate exp(x) near x=0 with [2/2] Padé approximant
    -- Power series: 1 + x + x^2/2 + x^3/6 + x^4/24
    let coeffs = [1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0]
    pade <- padeFromSeries coeffs 2 2 (-1.0) 1.0

    putStrLn "Approximation of exp(x):"
    let xs = [-0.5, -0.25, 0.0, 0.25, 0.5]
    ys <- evaluateMany pade xs
    let exacts = map exp xs
    let errors = zipWith (\y exact -> abs (y - exact)) ys exacts
    mapM_ (\(x, y, exact, err) ->
        printf "  x=%.2f: Padé=%.6f, exp(x)=%.6f, error=%.8f\n" x y exact err)
        (zip4 xs ys exacts errors)

    putStrLn ""

-- Example 5: With-style Usage
example5 :: IO ()
example5 = do
    putStrLn "5. With-style Usage (Bracket Pattern)"
    putStrLn (replicate 60 '-')

    -- Using withBernstein for automatic cleanup
    result <- withBernstein [1.0, 2.0, 6.0] 0.0 1.0 $ \poly -> do
        y <- evaluate poly 0.5
        dpoly <- derivative poly
        dy <- evaluate dpoly 0.5
        return (y, dy)

    let (y, dy) = result
    printf "f(0.5) = %.6f, f'(0.5) = %.6f\n" y dy

    -- Nested with usage
    result2 <- withBernstein [0.0, 1.0] 0.0 1.0 $ \num ->
               withBernstein [1.0, 2.0] 0.0 1.0 $ \den ->
               withRational num den $ \rat ->
                   evaluate rat 0.5

    printf "R(0.5) = %.6f\n" result2

    putStrLn ""
