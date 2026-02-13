module Test.Gelfgren (spec) where

import Test.Hspec
import Test.QuickCheck
import Control.Exception (evaluate, try, SomeException)
import Gelfgren

-- Tolerance for floating point comparisons
tolerance :: Double
tolerance = 1e-9

-- Helper to check approximate equality
(~=) :: Double -> Double -> Bool
x ~= y = abs (x - y) < tolerance

infix 4 ~=

spec :: Spec
spec = do
    describe "BernsteinPolynomial" $ do
        describe "creation" $ do
            it "creates a polynomial from coefficients" $ do
                poly <- bernsteinPolynomial [1.0, 2.0, 3.0] 0.0 1.0
                deg <- degree poly
                deg `shouldBe` 2

            it "throws error for empty coefficients" $ do
                result <- try $ bernsteinPolynomial [] 0.0 1.0
                result `shouldSatisfy` isLeft
              where
                isLeft :: Either SomeException a -> Bool
                isLeft (Left _) = True
                isLeft (Right _) = False

            it "throws error for invalid interval" $ do
                result <- try $ bernsteinPolynomial [1.0, 2.0] 1.0 0.0
                result `shouldSatisfy` isLeft
              where
                isLeft :: Either SomeException a -> Bool
                isLeft (Left _) = True
                isLeft (Right _) = False

        describe "evaluation" $ do
            it "evaluates constant polynomial correctly" $ do
                poly <- bernsteinPolynomial [2.0, 2.0] 0.0 1.0
                y <- evaluate poly 0.5
                y `shouldSatisfy` (~= 2.0)

            it "evaluates at left endpoint" $ do
                poly <- bernsteinPolynomial [1.0, 2.0, 3.0] 0.0 1.0
                y <- evaluate poly 0.0
                y `shouldSatisfy` (~= 1.0)

            it "evaluates at right endpoint" $ do
                poly <- bernsteinPolynomial [1.0, 2.0, 3.0] 0.0 1.0
                y <- evaluate poly 1.0
                y `shouldSatisfy` (~= 3.0)

            it "evaluates many points" $ do
                poly <- bernsteinPolynomial [1.0, 2.0, 3.0] 0.0 1.0
                ys <- evaluateMany poly [0.0, 0.5, 1.0]
                length ys `shouldBe` 3
                head ys `shouldSatisfy` (~= 1.0)
                last ys `shouldSatisfy` (~= 3.0)

        describe "derivative" $ do
            it "computes derivative of linear polynomial" $ do
                -- f(x) = 2x + 1, f'(x) = 2
                poly <- bernsteinPolynomial [1.0, 3.0] 0.0 1.0
                dpoly <- derivative poly
                dy <- evaluate dpoly 0.5
                dy `shouldSatisfy` (~= 2.0)

            it "derivative has correct degree" $ do
                poly <- bernsteinPolynomial [1.0, 2.0, 3.0] 0.0 1.0
                dpoly <- derivative poly
                deg <- degree dpoly
                deg `shouldBe` 1

        describe "integral" $ do
            it "computes integral of constant polynomial" $ do
                -- f(x) = 2, ∫f = 2x
                poly <- bernsteinPolynomial [2.0, 2.0] 0.0 1.0
                ipoly <- integral poly
                y0 <- evaluate ipoly 0.0
                y1 <- evaluate ipoly 1.0
                let area = y1 - y0
                area `shouldSatisfy` (~= 2.0)

            it "integral has correct degree" $ do
                poly <- bernsteinPolynomial [1.0, 2.0, 3.0] 0.0 1.0
                ipoly <- integral poly
                deg <- degree ipoly
                deg `shouldBe` 3

        describe "arithmetic" $ do
            it "adds two polynomials" $ do
                p1 <- bernsteinPolynomial [1.0, 1.0] 0.0 1.0
                p2 <- bernsteinPolynomial [2.0, 2.0] 0.0 1.0
                pSum <- add p1 p2
                y <- evaluate pSum 0.5
                y `shouldSatisfy` (~= 3.0)

            it "subtracts two polynomials" $ do
                p1 <- bernsteinPolynomial [3.0, 3.0] 0.0 1.0
                p2 <- bernsteinPolynomial [1.0, 1.0] 0.0 1.0
                pDiff <- subtract p1 p2
                y <- evaluate pDiff 0.5
                y `shouldSatisfy` (~= 2.0)

            it "multiplies two polynomials" $ do
                -- (x)(x) = x^2
                p1 <- bernsteinPolynomial [0.0, 1.0] 0.0 1.0
                p2 <- bernsteinPolynomial [0.0, 1.0] 0.0 1.0
                pProd <- multiply p1 p2
                deg <- degree pProd
                deg `shouldBe` 2
                y <- evaluate pProd 0.5
                y `shouldSatisfy` (~= 0.25)

            it "scales a polynomial" $ do
                poly <- bernsteinPolynomial [1.0, 2.0] 0.0 1.0
                scaled <- scale 3.0 poly
                y0 <- evaluate scaled 0.0
                y1 <- evaluate scaled 1.0
                y0 `shouldSatisfy` (~= 3.0)
                y1 `shouldSatisfy` (~= 6.0)

    describe "RationalFunction" $ do
        describe "creation" $ do
            it "creates a rational function" $ do
                num <- bernsteinPolynomial [0.0, 1.0] 0.0 1.0
                den <- bernsteinPolynomial [1.0, 1.0] 0.0 1.0
                rat <- rationalFunction num den
                y <- evaluate rat 0.5
                y `shouldSatisfy` (~= 0.5)

        describe "evaluation" $ do
            it "evaluates x/1 = x" $ do
                num <- bernsteinPolynomial [0.0, 1.0] 0.0 1.0
                den <- bernsteinPolynomial [1.0, 1.0] 0.0 1.0
                rat <- rationalFunction num den
                y <- evaluate rat 0.5
                y `shouldSatisfy` (~= 0.5)

            it "evaluates x/(1+x)" $ do
                num <- bernsteinPolynomial [0.0, 1.0] 0.0 1.0
                den <- bernsteinPolynomial [1.0, 2.0] 0.0 1.0
                rat <- rationalFunction num den
                y <- evaluate rat 0.5
                let expected = 0.5 / 1.5
                y `shouldSatisfy` (~= expected)

            it "evaluates at multiple points" $ do
                num <- bernsteinPolynomial [0.0, 1.0] 0.0 1.0
                den <- bernsteinPolynomial [1.0, 2.0] 0.0 1.0
                rat <- rationalFunction num den
                ys <- evaluateMany rat [0.0, 0.5, 1.0]
                length ys `shouldBe` 3
                head ys `shouldSatisfy` (~= 0.0)

    describe "PadeApproximant" $ do
        describe "creation" $ do
            it "creates a Padé approximant from series" $ do
                let coeffs = [1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0]
                pade <- padeFromSeries coeffs 2 2 (-1.0) 1.0
                y <- evaluate pade 0.0
                y `shouldSatisfy` (~= 1.0)

        describe "evaluation" $ do
            it "approximates exp(x) accurately" $ do
                let coeffs = [1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0]
                pade <- padeFromSeries coeffs 2 2 (-1.0) 1.0
                let x = 0.5
                y <- evaluate pade x
                let exact = exp x
                abs (y - exact) `shouldSatisfy` (< 1e-5)

            it "evaluates at multiple points" $ do
                let coeffs = [1.0, 1.0, 0.5, 1.0/6.0]
                pade <- padeFromSeries coeffs 2 1 (-1.0) 1.0
                ys <- evaluateMany pade [-0.5, 0.0, 0.5]
                length ys `shouldBe` 3
                ys !! 1 `shouldSatisfy` (~= 1.0)

    describe "With-style helpers" $ do
        it "withBernstein cleans up properly" $ do
            result <- withBernstein [1.0, 2.0] 0.0 1.0 $ \poly ->
                evaluate poly 0.5
            result `shouldSatisfy` (> 0)

        it "withRational works correctly" $ do
            num <- bernsteinPolynomial [0.0, 1.0] 0.0 1.0
            den <- bernsteinPolynomial [1.0, 1.0] 0.0 1.0
            result <- withRational num den $ \rat ->
                evaluate rat 0.5
            result `shouldSatisfy` (~= 0.5)

        it "withPade works correctly" $ do
            result <- withPade [1.0, 1.0, 0.5] 1 1 (-1.0) 1.0 $ \pade ->
                evaluate pade 0.0
            result `shouldSatisfy` (~= 1.0)

    describe "Error handling" $ do
        it "catches null pointer errors" $ do
            result <- try $ bernsteinPolynomial [] 0.0 1.0
            result `shouldSatisfy` isLeft
          where
            isLeft :: Either SomeException a -> Bool
            isLeft (Left _) = True
            isLeft (Right _) = False
