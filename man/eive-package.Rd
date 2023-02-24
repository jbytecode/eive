\name{eive-package}
\alias{eive-package}
\alias{eive}
\docType{package}
\title{
EIVE - Errors-in-Variable estimation
}
\description{
This package includes functions for compact genetic algorithms and errors-in-variable estimation. The function 'eive' performs a genetic search to reduce the errors-in-variable bias in ordinary least squares estimator. 

Change log:

# In version 3.1.2

- Add eive.cga.formula for lm() compatible regression settings.

- Implement roxygen type documentation system

# In version 3.1.1, 

- Enhance document for multivariate eive (eivem)

# In version 3.1, 

- Errors-in-variables with multiple response variables.

- eive.cga() now returns 2 new items: $cleanx and $measurementerror which are also accessable using the keys $proxy$fitted.values and $proxy$residuals, respectively.

- eivem() added for multiple y values

- new tests added and can be triggered using devtools:::test()

- updated docs

# In version 2.1, more speed improvements by using lm.fit instead lm in critical code. 

# In version 2.0, some routines are rewritten in C++ and wrapped using Rcpp so a substantial speed improvement achieved.

}
\details{
\tabular{ll}{
Package: \tab eive\cr
Type: \tab Package\cr
Version: \tab 2.1\cr
Date: \tab 2014-07-04\cr
License: \tab GPL\cr
}
}
\author{
Mehmet Hakan Satman <mhsatman@istanbul.edu.tr>
Erkin Diyarbakirlioglu <ediyarbakirlioglu@gmail.com>


Maintainer: Mehmet Hakan Satman <mhsatman@istanbul.edu.tr>
}
