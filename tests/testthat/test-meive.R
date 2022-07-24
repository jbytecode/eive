library(testthat)
library(eive)

test_that("Multivariate Eive - three y variables - n = 30", {
    set.seed(12345)
    n <- 30
    clean_x1 <- rnorm(n, mean = 10, sd = sqrt(7))
    delta_x1 <- rnorm(n, mean = 0, sd = sqrt(3))

    dirty_x1 <- clean_x1 + delta_x1

    real_beta_0 <- 20
    real_beta_1 <- 10
    real_coefs <- c(real_beta_0, real_beta_1)

    e1 <- rnorm(n, mean = 0, sd = sqrt(5))
    e2 <- rnorm(n, mean = 0, sd = sqrt(5))
    e3 <- rnorm(n, mean = 0, sd = sqrt(5))

    y1 <- real_beta_0 + real_beta_1 * clean_x1 + e1
    y2 <- real_beta_0 + real_beta_1 * clean_x1 + e2
    y3 <- real_beta_0 + real_beta_1 * clean_x1 + e3

    result <- eivem(
        dirtyx = dirty_x1,
        otherx = NULL,
        numdummies = 10,
        popsize = 20,
        y = cbind(y1, y2, y3)
    )

    ols <- result$ols[[1]]
    eive <- result$eive[[1]]

    ols_coefs <- ols$coefficients
    eive_coefs <- eive$coefficients

    distance1 <- sum((ols_coefs - real_coefs)^2)
    distance2 <- sum((eive_coefs - real_coefs)^2)

    expect_true(
        distance1 > distance2
    )

    result_names <- names(result)
    expect_true("ols" %in% result_names)
    expect_true("eive" %in% result_names)
    expect_true("proxy" %in% result_names)
    expect_true("cleanedx" %in% result_names)
    expect_true("measurementerror" %in% result_names)
})