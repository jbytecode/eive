library("testthat")

test_that("Multiple Eive - two variables - n = 30", {
    set.seed(12345)
    n <- 30
    clean_x1 <- rnorm(n, mean = 10, sd = sqrt(7))
    clean_x2 <- rnorm(n, mean = 10, sd = sqrt(7))
    delta_x1 <- rnorm(n, mean = 0, sd = sqrt(3))
    delta_x2 <- rnorm(n, mean = 0, sd = sqrt(3))
    e <- rnorm(n, mean = 0, sd = sqrt(5))

    real_beta_0 <- 20
    real_beta_1 <- 10
    real_beta_2 <- 10
    real_coefs <- c(real_beta_0, real_beta_1, real_beta_2)

    y <- real_beta_0 + real_beta_1 * clean_x1 + real_beta_2 * clean_x2 + e

    dirty_x1 <- clean_x1 + delta_x1
    dirty_x2 <- clean_x2 + delta_x2

    result <- meive.cga(
        dirtyxmat = cbind(dirty_x1, dirty_x2),
        otherx = NULL,
        numdummiesforeach = 10,
        popsize = 60,
        y = y
    )

    ols <- result$ols
    eive <- result$eive

    ols_coefs <- ols$coefficients
    eive_coefs <- eive$coefficients

    # print("OLS:")
    # print(ols_coefs)
    # print("MEIVE:")
    # print(eive_coefs)

    distance1 <- sum((ols_coefs - real_coefs)^2)
    distance2 <- sum((eive_coefs - real_coefs)^2)

    expect_true(
        distance1 > distance2
    )
})