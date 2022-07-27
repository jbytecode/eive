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
    probs <- rep(0.5, 1000)
    v <- rep(0, 1000)
    result <- cga_generate_chromosome(probs, v)
    
    # Function returns NULL
    expect_null(result)
    
    # Not all of the values are same
    expect_false(all(v == 0))
    expect_false(all(v == 1))
})