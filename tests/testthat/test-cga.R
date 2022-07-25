library("testthat")


test_that("CGA - All 10 bits are 1s", {
    f <- function(bits) {
        return(-sum(bits))
    }
    n <- 10
    result <- cga(n, 100, f)
    expect_equal(result, rep(1, n))
})

test_that("CGA - All 100 bits are 1s", {
    f <- function(bits) {
        return(-sum(bits))
    }
    n <- 100
    result <- cga(n, 100, f)
    expect_equal(result, rep(1, n))
})

test_that("CGA - All 100 bits are 0s", {
    f <- function(bits) {
        return(sum(bits))
    }
    n <- 100
    result <- cga(n, 100, f)
    expect_equal(result, rep(0, n))
})

test_that("CGA - Generate Chromosome", {
    probs <- rep(0.5, 10)
    v <- rep(0, 10)
    cga_generate_chromosome(probs, v)
    for (element in v) {
        expect_gte(element, 0.0)
        expect_lte(element, 1.0)
    }
})