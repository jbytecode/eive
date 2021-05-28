# Generates one or two exploratory variables linear model
# with first one is subject to error
generate.eive.data <- function(n,
                               e.sd,
                               delta.sd,
                               seed = 12345,
                               useotherx = FALSE) {
    set.seed(seed)
    e <- rnorm(n, 0, e.sd)
    x <- rnorm(n)
    x1 <- NULL
    y <- NULL
    if (useotherx) {
        x1 <- rnorm(n)
    }
    delta <- rnorm(n, 0, delta.sd)
    xdelta <- x + delta
    if (useotherx) {
        y <- 5 + 5 * x + 5 * x1 + e
    } else {
        y <- 5 + 5 * x + e
    }
    data <- cbind(xdelta, x1, y)
    return(data)
}



eive.cga <- function(dirtyx,
                     otherx = NULL,
                     y,
                     numdummies = 10,
                     popsize = 20) {
    ols_dirty <- NULL
    ols_proxy <- NULL
    ols_best <- NULL
    n <- length(y)
    f <- function(d) {
        ols <- NULL
        m <- matrix(d, nrow = n)
        ols_proxy <- lm.fit(cbind(1, m), dirtyx)
        x_proxy <- ols_proxy$fitted.values
        if (is.null(otherx)) {
            ols <- lm.fit(cbind(1, x_proxy), y)
        } else {
            ols <- lm.fit(cbind(1, x_proxy, otherx), y)
        }
        return(sum(ols$residuals^2))
    }

    ga <- cga(
        evalFunc = f,
        chsize = n * numdummies,
        popsize = popsize
    )
    # best<-cga_generate_chromosome(ga)
    best <- as.integer(ga)
    if (is.null(otherx)) {
        ols_dirty <- lm(y ~ dirtyx)
    } else {
        ols_dirty <- lm(y ~ dirtyx + otherx)
    }

    ols_proxy <- lm(dirtyx ~ matrix(best, nrow = n))

    if (is.null(otherx)) {
        ols_best <- lm(y ~ ols.proxy$fitted.values)
    } else {
        ols_best <- lm(y ~ ols.proxy$fitted.values + otherx)
    }

    result <- list(ols = ols_dirty, eive = ols_best, proxy = ols_proxy)
    return(result)
}



meive.cga <- function(dirtyxmat,
                      otherx = NULL,
                      y,
                      numdummiesforeach = 10,
                      popsize = 20) {
    ols_dirty <- NULL
    ols_best <- NULL
    num_dirtyx <- ncol(dirtyxmat)
    num_dummies <- numdummiesforeach * num_dirtyx
    n <- length(y)

    get_x_proxies <- function(bits) {
        bits_mat <- matrix(bits, nrow = n)
        m <- 1
        x_proxies <- matrix(0, ncol = num_dirtyx, nrow = n)
        for (i in 1:num_dirtyx) {
            bits_submat <- bits_mat[, m:(m + numdummiesforeach - 1)]
            ols_proxy <- lm.fit(cbind(1, bits_submat), dirtyxmat[, i])
            x_proxy <- ols_proxy$fitted.values
            x_proxies[, i] <- x_proxy
            m <- m + numdummiesforeach
        }
        return(x_proxies)
    }

    cost <- function(bits) {
        ols <- NULL
        x_proxies <- get_x_proxies(bits)
        if (is.null(otherx)) {
            ols <- lm.fit(cbind(1, x_proxies), y)
        } else {
            ols <- lm.fit(cbind(1, x_proxies, otherx), y)
        }
        return(sum(ols$residuals^2))
    }

    ga <- cga(
        evalFunc = cost,
        chsize = n * numdummiesforeach * num_dirtyx,
        popsize = popsize
    )

    best_bits <- as.integer(ga)

    if (is.null(otherx)) {
        ols_dirty <- lm(y ~ dirtyxmat)
    } else {
        ols_dirty <- lm(y ~ cbind(dirtyxmat, otherx))
    }

    clean_x_mat <- get_x_proxies(best_bits)
    if (is.null(otherx)) {
        ols_best <- lm(y ~ clean_x_mat)
    } else {
        ols_best <- lm(y ~ cbind(clean_x_mat, otherx))
    }

    result <- list(
        ols = ols_dirty,
        eive = ols_best,
        corrected_x_mat = clean_x_mat
    )

    return(result)
}