library(testthat)
library(eive)

test_that("eive.cga classical example using formula", {
    tol <- 0.01

    euclidean <- function(u, v) {
        return(sqrt(sum((u - v)^2)))
    }

    require("eive")

    set.seed(12345)

    n <- 30

    clean_x <- rnorm(n, mean = 10, sd = sqrt(7))

    delta_x <- rnorm(n, mean = 0, sd = sqrt(3))

    e <- rnorm(n, mean = 0, sd = sqrt(5))

    y <- 20 + 10 * clean_x + e

    dirty_x <- clean_x + delta_x

    mydata <- data.frame(y = y, dirtyx = dirty_x)

    result <- eive.cga.formula(
        formula = y ~ dirtyx,
        dirtyx.varname = "dirtyx",
        data = mydata,
        numdummies = 10
    )

    par_real <- c(20, 10)
    par_ols <- result$ols$coefficients
    par_eive <- result$eive$coefficients

    dist1 <- euclidean(par_real, par_ols)
    dist2 <- euclidean(par_real, par_eive)

    expect_true(dist2 < dist1)

    expect_equal(as.vector(par_eive), c(23.863, 9.229), tolerance = tol)
    expect_equal(as.vector(par_ols), c(63.590, 5.533), tolerance = tol)

    L <- length(result$measurementerror)
    for (i in 1:L) {
        expect_equal(result$proxy$residuals[i], result$measurementerror[i])
    }

    result_names <- names(result)
    expect_true("ols" %in% result_names)
    expect_true("eive" %in% result_names)
    expect_true("proxy" %in% result_names)
    expect_true("cleanedx" %in% result_names)
    expect_true("measurementerror" %in% result_names)
})




