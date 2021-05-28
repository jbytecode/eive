library("testthat")

test_status <- function(str_message) {
    cat("* Doing test: ", str_message, "\n")
}

test_that("CGA - All 10 bits are 1s", {
    test_status("CGA - All 10 bits are 1s")
    f <- function(bits) {
        return(-sum(bits))
    }
    n <- 10
    result <- cga(n, 100, f)
    expect_equal(result, rep(1, n))
})

test_that("CGA - All 100 bits are 1s", {
    test_status("CGA - All 100 bits are 1s")
    f <- function(bits) {
        return(-sum(bits))
    }
    n <- 100
    result <- cga(n, 100, f)
    expect_equal(result, rep(1, n))
})

test_that("CGA - All 100 bits are 0s", {
    test_status("CGA - All 100 bits are 0s")
    f <- function(bits) {
        return(sum(bits))
    }
    n <- 100
    result <- cga(n, 100, f)
    expect_equal(result, rep(0, n))
})