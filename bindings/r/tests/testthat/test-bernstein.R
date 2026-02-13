test_that("BernsteinPolynomial creation works", {
    p <- bernstein_polynomial(c(1, 2, 3), 0, 1)
    expect_s3_class(p, "BernsteinPolynomial")
    expect_equal(p$degree(), 2)
})

test_that("BernsteinPolynomial evaluation works", {
    p <- bernstein_polynomial(c(1, 2, 3), 0, 1)

    # Evaluate at endpoints
    expect_equal(p$evaluate(0.0)[1], 1.0)
    expect_equal(p$evaluate(1.0)[1], 3.0)

    # Evaluate at multiple points
    x <- c(0, 0.5, 1)
    y <- p$evaluate(x)
    expect_length(y, 3)
    expect_true(all(is.finite(y)))
})

test_that("Polynomial arithmetic works", {
    p <- bernstein_polynomial(c(1, 2), 0, 1)
    q <- bernstein_polynomial(c(3, 4), 0, 1)

    sum <- p$add(q)
    expect_s3_class(sum, "BernsteinPolynomial")

    diff <- p$subtract(q)
    expect_s3_class(diff, "BernsteinPolynomial")

    prod <- p$multiply(q)
    expect_s3_class(prod, "BernsteinPolynomial")
})

test_that("Derivative computation works", {
    p <- bernstein_polynomial(c(1, 2, 3), 0, 1)
    p_prime <- p$derivative()

    expect_s3_class(p_prime, "BernsteinPolynomial")
    expect_equal(p_prime$degree(), 1)
})

test_that("Degree elevation works", {
    p <- bernstein_polynomial(c(1, 2), 0, 1)
    elevated <- p$elevate()

    expect_equal(elevated$degree(), p$degree() + 1)

    # Values should be preserved
    x <- 0.5
    expect_equal(p$evaluate(x), elevated$evaluate(x), tolerance = 1e-10)
})

test_that("Invalid inputs throw errors", {
    expect_error(bernstein_polynomial(numeric(0), 0, 1))
    expect_error(bernstein_polynomial(c(1, 2), 1, 0))  # a > b
})
