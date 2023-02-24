library(compiler)

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


eive.cga.formula <- function(
    formula,
    data, 
    dirtyx.varname,
    numdummies = 10,
    popsize = 20
){
    if (!is.data.frame(data)) {
        stop("data should be in type of data.frame")
    }
    if (!is.language(formula)) {
        stop("formula should be in form of var ~ var + var ...")
    }
    if (!is.character(dirtyx.varname)) {
        stop("dirtyx.varname should be in type of string")
    }
    designmat <- as.data.frame(model.matrix(object = formula, data = data))
    responsevar <- model.frame(formula = formula, data = data)
    y <- responsevar[, 1]
    n <- length(y)
    p <- ncol(designmat)
    dirtyx <- as.vector(data[, dirtyx.varname])
    otherx <- designmat[, setdiff(names(designmat), dirtyx.varname)]

    if (is.vector(otherx)){
        otherx <- as.matrix(otherx)
    }
    has.intercept <- all.equal(otherx[, 1], rep(1, n))
    if (has.intercept) {
        if (ncol(otherx) > 1) {
            otherx <- as.matrix(otherx[, 2:(ncol(otherx))])
        }else{
            otherx <- NULL
        }
    }

    return(
        eive.cga(
            dirtyx = dirtyx, 
            otherx = otherx,
            y = y, 
            numdummies = numdummies, 
            popsize = popsize
        ))
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
    
    best <- as.integer(ga)

    if (is.null(otherx)) {
        ols_dirty <- lm(y ~ dirtyx)
    } else {
        ols_dirty <- lm(y ~ dirtyx + otherx)
    }

    ols_proxy <- lm(dirtyx ~ matrix(best, nrow = n))

    if (is.null(otherx)) {
        ols_best <- lm(y ~ ols_proxy$fitted.values)
    } else {
        ols_best <- lm(y ~ ols_proxy$fitted.values + otherx)
    }

    cleanedx <- ols_proxy$fitted.values
    measurementerror <- ols_proxy$residuals
    
    result <- list(
      ols = ols_dirty,
      eive = ols_best,
      proxy = ols_proxy,
      cleanedx = cleanedx,
      measurementerror = measurementerror
    )
    
    return(result)
}


