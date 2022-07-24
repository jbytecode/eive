\name{eivem}
\alias{eivem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Errors-in-variable estimation in linear regression with compact genetic algorithms
}
\usage{
eivem(dirtyx, otherx = NULL, y, numdummies = 10, popsize = 20)
}
\arguments{
  \item{dirtyx}{
Vector of values of independent variable measured with error
}
  \item{otherx}{
Matrix of other independent variables.
}
  \item{y}{
Matrix of values of dependent variables. In column represent a single Y, each row represents an observation. So y is an nxp matrix of p dependent variables.
}
  \item{numdummies}{
Number of dummy variables used in algorithm. By default, it is 10.
}
  \item{popsize}{
Population size parameter used in CGA. By default it is 20.
}
}
\details{
Algorithm performs a genetic search to separate mismeasured independent variable into clean and error parts.
}
\value{
 \item{ols}{lm object calculated using original values}
 \item{eive}{lm object calculated using the predicted variable by eive}
 \item{proxy}{lm object of proxy regression obtained by genetic search.}
 \item{cleanedx}{Error-free estimate of the x variable (dirtyx) that is measured with error.}
 \item{measurementerror}{Estimate of the measurement error.}
}
\author{
Mehmet Hakan Satman <mhsatman@istanbul.edu.tr>
Erkin Diyarbakirlioglu <ediyarbakirlioglu@gmail.com>
}
\description{
Method performs a genetic search to find dummy variables that used in a two stage linear regression to reduce errors-in-variables
bias in linear regression. This is the multivariate version of the CGA based EIVE algorithm.
}

\examples{
# Creating an artificial data

# Loading required package
require ("eive")

# Setting random number generator seed to 12345
# so each time the script runs, same numbers will
# be generated
set.seed(12345)

# Number of observations is set to 30
n <- 30

# Unobserved X values are drawn from a Normal distribution
# with mean 10 and variance 7
clean_x1 <- rnorm(n, mean=10, sd=sqrt(7))
clean_x2 <- rnorm(n, mean=10, sd=sqrt(7))

# Measurement error values are dranw from a Normal distribution
# with mean 0 and variance 3
delta_x1 <- rnorm(n, mean=0,  sd=sqrt(3))

# Error term of regression. Normally distributed with mean 0 and 
# variance 5
e <- rnorm(n, mean=0, sd=sqrt(5))

# Generating Y values using the linear model
# In this model, intercept is 20 and slope is 10. 
y1 <- 20 + 10 * clean_x1 + 10 * clean_x2 + e
y2 <- 10 + 5 * clean_x1  + 5 * clean_x2 + e

# Generating observed X values by adding measurement errors
# to unobserved X
dirty_x1 <- clean_x1 + delta_x1

# Performs a genetic search to find dummy variables that
# used in two stage least squares.
# Please un-comment the line below
result <- eivem (dirtyx = dirty_x1, y=cbind(y1, y2), numdummies=10)

# Print the result
# Please un-comment the line below
# print(result)

########################################### OUTPUT #############################################
#> result
#$ols
#$ols[[1]]
#
#Call:
#lm(formula = y[, reg.index] ~ dirtyx)
#
#Coefficients:
#(Intercept)       dirtyx  
#    156.975        6.696  
#
#
#$ols[[2]]
#
#Call:
#lm(formula = y[, reg.index] ~ dirtyx)
#
#Coefficients:
#(Intercept)       dirtyx  
#     79.826        3.229  
#
#
#
#$eive
#$eive[[1]]
#
#Call:
#lm(formula = y[, reg.index] ~ ols_proxy$fitted.values)
#
#Coefficients:
#            (Intercept)  ols_proxy$fitted.values  
#                  20.05                    19.19  
#
#
#$eive[[2]]
#
#Call:
#lm(formula = y[, reg.index] ~ ols_proxy$fitted.values)
#
#Coefficients:
#            (Intercept)  ols_proxy$fitted.values  
#                 11.451                    9.466  
#
#
#
#$proxy
#
#Call:
#lm(formula = dirtyx ~ matrix(best, nrow = n))
#
#Coefficients:
#             (Intercept)   matrix(best, nrow = n)1   matrix(best, nrow = n)2  
#                 9.60535                   0.51227                   3.35246  
# matrix(best, nrow = n)3   matrix(best, nrow = n)4   matrix(best, nrow = n)5  
#                 0.05283                  -0.87955                  -0.62274  
# matrix(best, nrow = n)6   matrix(best, nrow = n)7   matrix(best, nrow = n)8  
#                -0.52760                   1.50575                  -0.68810  
# matrix(best, nrow = n)9  matrix(best, nrow = n)10  
#                 0.27033                   0.37630  
#
#
#
######################################### END OF OUTPUT ##########################################

}
