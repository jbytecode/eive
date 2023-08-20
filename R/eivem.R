library(compiler)


#' Performs CGA based errors-in-variables correction for a given set of variables
#' in case of multiple Y variables are provided.
#' @description  A single independent variable is supposed to be measured subject to error.
#'               This functions is the multivariate version of the classical algorithm.
#'               Additional response variables are used to get better estimates.
#'
#' @param dirtyx Vector of independent variable that is measured with error.
#' @param otherx Matrix of other independent variables. If the model has
#'        a single independent variable, it is NULL by default.
#' @param y Matrix of response variables. Y_i is placed in the ith row of the matrix.
#' @param numdummies Number of dummy variables used in auxiliary regression. Default is 10.
#' @param popsize Population size parameter for compact genetic algorithm. Default is 20.
#' 1/popsize is the mutation rate.
#' @return A list() of regression equations.
#' @slot ols List of lm objects calculated using original values
#' @slot eive List of lm objects calculated using the predicted variable by eive
#' @slot proxy lm object of proxy regression obtained by genetic search.
#' @slot cleanedx Error-free estimate of the x variable (dirtyx)
#'       that is measured with error.
#' @slot measurementerror Estimate of the measurement error.
#' @export 
#' @examples
#' # Creating an artificial data
#' 
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
#' clean_x1 <- rnorm(n, mean = 10, sd = sqrt(7))
#' clean_x2 <- rnorm(n, mean = 10, sd = sqrt(7))
#' 
#' # Measurement error values are dranw from a Normal distribution
#' # with mean 0 and variance 3
#' delta_x1 <- rnorm(n, mean = 0, sd = sqrt(3))
#' 
#' # Error term of regression. Normally distributed with mean 0 and
#' # variance 5
#' e1 <- rnorm(n, mean = 0, sd = sqrt(5))
#' e2 <- rnorm(n, mean = 0, sd = sqrt(5))
#' 
#' # Generating Y values using the linear model
#' # In this model, intercept is 20 and slope is 10.
#' y1 <- 20 + 10 * clean_x1 + 10 * clean_x2 + e1
#' y2 <- 10 + 5 * clean_x1 + 5 * clean_x2 + e2
#' 
#' # Generating observed X values by adding measurement errors
#' # to unobserved X
#' dirty_x1 <- clean_x1 + delta_x1
#' 
#' # Performs a genetic search to find dummy variables that
#' # used in two stage least squares.
#' # Please un-comment the line below
#' result <- eivem(dirtyx = dirty_x1, otherx = clean_x2, y = cbind(y1, y2), numdummies = 10)
#' 
#' # Print the result
#' # Please un-comment the line below
#' # print(result)
#' 
#' ########################################### OUTPUT #############################################
#' #> result
#' # $ols
#' # $ols[[1]]
#' #
#' # Call:
#' # lm(formula = y[, reg.index] ~ dirtyx + otherx)
#' #
#' # Coefficients:
#' # (Intercept)       dirtyx       otherx
#' #      54.141        6.067       10.137
#' #
#' #
#' # $ols[[2]]
#' #
#' # Call:
#' # lm(formula = y[, reg.index] ~ dirtyx + otherx)
#' #
#' # Coefficients:
#' # (Intercept)       dirtyx       otherx
#' #      24.814        3.205        5.089
#' #
#' #
#' #
#' # $eive
#' # $eive[[1]]
#' #
#' # Call:
#' # lm(formula = y[, reg.index] ~ ols_proxy$fitted.values + otherx)
#' #
#' # Coefficients:
#' #             (Intercept)  ols_proxy$fitted.values                   otherx
#' #                  24.737                    9.727                    9.147
#' #
#' #
#' # $eive[[2]]
#' #
#' # Call:
#' # lm(formula = y[, reg.index] ~ ols_proxy$fitted.values + otherx)
#' #
#' # Coefficients:
#' #             (Intercept)  ols_proxy$fitted.values                   otherx
#' #                   8.313                    5.240                    4.552
#' #
#' #
#' 
#' # $proxy
#' #
#' # Call:
#' # lm(formula = dirtyx ~ matrix(best, nrow = n))
#' #
#' # Coefficients:
#' #              (Intercept)   matrix(best, nrow = n)1   matrix(best, nrow = n)2
#' #                 6.314397                 -0.211580                  1.729143
#' #  matrix(best, nrow = n)3   matrix(best, nrow = n)4   matrix(best, nrow = n)5
#' #                 1.994915                  0.947531                 -0.363107
#' #  matrix(best, nrow = n)6   matrix(best, nrow = n)7   matrix(best, nrow = n)8
#' #                 0.001768                  1.742553                 -0.023750
#' #  matrix(best, nrow = n)9  matrix(best, nrow = n)10
#' #                 0.134750                  2.324853
#' #
#' #
#' # $cleanedx
#' #         1         2         3         4         5         6         7         8
#' # 12.730307 12.130102 11.065586  9.795474 12.697138  6.450915 12.673388 10.516553
#' #         9        10        11        12        13        14        15        16
#' # 11.095771  7.981887 11.694464 14.841812 11.098755 12.290371  8.988344 12.704789
#' #        17        18        19        20        21        22        23        24
#' #  7.671861  9.477178 13.458999 10.964004 11.465852 14.591473  9.771724  6.239335
#' #        25        26        27        28        29        30
#' #  6.425397 15.031410  8.992839 12.808138 13.435249  9.799758
#' #
#' # $measurementerror
#' #          1          2          3          4          5          6          7
#' # -0.9220426 -2.5783644 -0.3964263  1.7585818 -2.1106159 -4.4345451  0.5319987
#' #          8          9         10         11         12         13         14
#' #  1.5127360 -0.9523682 -2.6583539 -1.9074299 -1.3927085 -1.9356982  3.1225578
#' #         15         16         17         18         19         20         21
#' #  1.4554922  1.0891572  1.4141792 -1.7600789  0.3310142  1.5952156  1.7146703
#' #         22         23         24         25         26         27         28
#' #  1.0669497 -2.0036393  3.9419318  1.0296643  2.9783401  0.8968531 -1.7001587
#' #         29         30
#' # -0.8864360  1.1995241
#' #
#' ######################################### END OF OUTPUT ##########################################
eivem <- function(dirtyx, otherx = NULL, y, numdummies = 10, popsize = 20) {
  ols_dirty <- list()
  ols_proxy <- NULL
  ols_best <- list()
  my_dim <- dim(y)
  n <- my_dim[1]
  yp <- my_dim[2]

  f <- function(d) {
    ols <- NULL
    m <- matrix(d, nrow = n)
    ols_proxy <- lm.fit(cbind(1, m), dirtyx)
    x_proxy <- ols_proxy$fitted.values
    square_sums <- 0
    for (reg.index in 1:yp) {
      if (is.null(otherx)) {
        ols <- lm.fit(cbind(1, x_proxy), y[, reg.index])
      } else {
        ols <- lm.fit(cbind(1, x_proxy, otherx), y[, reg.index])
      }
      square_sums <- square_sums + sum(ols$residuals^2)
    }
    return(square_sums)
  }

  g <- cmpfun(f)

  ga <- cga(evalFunc = g, chsize = n * numdummies, popsize = popsize)

  best <- as.integer(ga)

  for (reg.index in 1:yp) {
    if (is.null(otherx)) {
      ols_dirty[[reg.index]] <- lm(y[, reg.index] ~ dirtyx)
    } else {
      ols_dirty[[reg.index]] <- lm(y[, reg.index] ~ dirtyx + otherx)
    }
    ols_proxy <- lm(dirtyx ~ matrix(best, nrow = n))
    if (is.null(otherx)) {
      ols_best[[reg.index]] <- lm(y[, reg.index] ~ ols_proxy$fitted.values)
    } else {
      ols_best[[reg.index]] <- lm(
        y[, reg.index] ~ ols_proxy$fitted.values + otherx
      )
    }
  } # end of for loop

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
# End of function
