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


#' Performs CGA based errors-in-variables correction for given formula and data.
#' A single independent variable is supposed to be measured subject to error.
#'
#' @param formula Formula object.
#' @param data data.frame that holds the regression data.
#' @param dirtyx.varname String key value of the erroneous independent variable.
#' @param numdummies Number of dummy variables used in auxiliary regression.
#' @param popsize Population size parameter for compact genetic algorithm. 
#' 1/popsize is the mutation rate.
#' @return A list() of regression equations.
#' @slot ols lm object calculated using original values
#' @slot eive lm object calculated using the predicted variable by eive
#' @slot proxy lm object of proxy regression obtained by genetic search.
#' @slot cleanedx Error-free estimate of the x variable (dirtyx) 
#'       that is measured with error.
#' @slot measurementerror Estimate of the measurement error.
#' @export
#' @examples
#' set.seed(12345)
#' n <- 30
#' clean_x <- rnorm(n, mean = 10, sd = sqrt(7))
#' delta_x <- rnorm(n, mean = 0, sd = sqrt(3))
#' 
#' e <- rnorm(n, mean = 0, sd = sqrt(5))
#' y <- 20 + 10 * clean_x + e
#' 
#' dirty_x <- clean_x + delta_x
#' 
#' mydata <- data.frame(y = y, dirtyx = dirty_x)
#' 
#' result <- eive.cga.formula(
#'      formula = y ~ dirtyx,
#'      dirtyx.varname = "dirtyx",
#'      data = mydata,
#'      numdummies = 10
#' )
#' @seealso eive.cga
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





#' Performs CGA based errors-in-variables correction for a given set of variables.
#' A single independent variable is supposed to be measured subject to error.
#'
#' @param dirtyx Vector of independent variable that is measured with error.
#' @param otherx Matrix of other independent variables. If the model has 
#'        a single independent variable, it is NULL by default.
#' @param y Vector of response variable
#' @param numdummies Number of dummy variables used in auxiliary regression.
#' @param popsize Population size parameter for compact genetic algorithm.
#' 1/popsize is the mutation rate.
#' @return A list() of regression equations.
#' @slot ols lm object calculated using original values
#' @slot eive lm object calculated using the predicted variable by eive
#' @slot proxy lm object of proxy regression obtained by genetic search.
#' @slot cleanedx Error-free estimate of the x variable (dirtyx)
#'       that is measured with error.
#' @slot measurementerror Estimate of the measurement error.
#' @export 
#' @examples 
#' # Creating an artificial data
#' # Loading required package
#' require("eive")
#' 
#' # Setting random number generator seed to 12345
#' # so each time the script runs, same numbers will
#' # be generated
#' set.seed(12345)
#' 
#' # Number of observations is set to 30
#' n <- 30
#' 
#' # Unobserved X values are drawn from a Normal distribution
#' # with mean 10 and variance 7
#' clean.x <- rnorm(n, mean = 10, sd = sqrt(7))
#' 
#' # Measurement error values are dranw from a Normal distribution
#' # with mean 0 and variance 3
#' delta.x <- rnorm(n, mean = 0, sd = sqrt(3))
#' 
#' # Error term of regression. Normally distributed with mean 0 and
#' # variance 5
#' e <- rnorm(n, mean = 0, sd = sqrt(5))
#' 
#' # Generating Y values using the linear model
#' # In this model, intercept is 20 and slope is 10.
#' y <- 20 + 10 * clean.x + e
#' 
#' # Generating observed X values by adding measurement errors
#' # to unobserved X
#' dirty.x <- clean.x + delta.x
#' 
#' # Performs a genetic search to find dummy variables that
#' # used in two stage least squares.
#' # Please un-comment the line below
#' # result <- eive.cga (dirtyx=dirty.x, y=y, numdummies=10)
#' 
#' # Print the result
#' # Please un-comment the line below
#' # print(result)
#' 
#' ########################################### OUTPUT #############################################
#' # $ols
#' #
#' # Call:
#' # lm(formula = y ~ dirtyx)
#' #
#' # Coefficients:
#' # (Intercept)       dirtyx
#' #     63.590        5.533
#' #
#' #
#' # $eive
#' #
#' # Call:
#' # lm(formula = y ~ ols.proxy$fitted.values)
#' #
#' # Coefficients:
#' #            (Intercept)  ols.proxy$fitted.values
#' #                 23.863                    9.229
#' #
#' #
#' # $proxy
#' #
#' # Call:
#' # lm(formula = dirtyx ~ matrix(best, nrow = n))
#' #
#' # Coefficients:
#' #              (Intercept)   matrix(best, nrow = n)1   matrix(best, nrow = n)2
#' #                 12.9321                   -0.6252                   -1.9923
#' # matrix(best, nrow = n)3   matrix(best, nrow = n)4   matrix(best, nrow = n)5
#' #                  0.7537                   -0.7076                   -0.5247
#' # matrix(best, nrow = n)6   matrix(best, nrow = n)7   matrix(best, nrow = n)8
#' #                 -0.9196                   -2.0802                   -0.9246
#' # matrix(best, nrow = n)9  matrix(best, nrow = n)10
#' #                 -0.6164                    1.9694
#' ######################################### END OF OUTPUT ##########################################
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


