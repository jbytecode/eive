# eive
An R package for Errors-in-variables estimation in linear regression

## Installation

### Install stable version from CRAN

```R
install.packages("eive")
```


### Install development version 

Please install ```devtools``` package before installing ```eive```:

```R
install.packages("devtools")
```

then install the package from the github repo using

```R
devtools::install_github(repo = "https://github.com/jbytecode/eive") 
```

# The Problem 

Suppose the linear regression model is 

$$
y = \beta_0 + \beta_1 x^* + \varepsilon
$$

where $y$ is n-vector of the response variable, $\beta_0$ and $\beta_1$ are unknown regression parameteres, $\varepsilon$ is the iid. error term, $x^*$ is the unknown n-vector of the independent variable, and $n$ is the number of observations.

We call $x^*$ unknown because in some situations the true values of the variable cannot be visible or directly observable, or observable with some measurement error. Now suppose that $x$ is the observable version of the true values and it is defined as 

$$
x = x^* + \delta
$$

where $\delta$ is the measurement error and $x$ is the erroneous version of the true $x^*$. If the estimated model is 

$$
\hat{y} = \hat{\beta_0} + \hat{\beta_1}x 
$$

then the ordinary least squares (OLS) estimates are no longer unbiased and even consistent. 

Eive-cga is an estimator devised for this problem. The aim is to reduce the errors-in-variable bias with some cost of increasing the variance. At the end, the estimator obtains lower Mean Square Error (MSE) values defined as

$$
MSE(\hat{\beta_1}) = Var(\hat{\beta_1}) + Bias^2(\hat{\beta_1})
$$

for the Eive-cga estimator. For more detailed comparisons, see the original paper given in the Citation part. 

# Usage 

For the single variable case 

```R 
> eive(dirtyx = dirtyx, y = y, otherx = nothing) 
```

and for the multiple regression 

```R 
> eive(dirtyx = dirtyx, y = y, otherx = matrixofotherx) 
```

and for the multiple regression with formula object 

```R 
> eive(formula = y ~ x1 + x2 + x3, dirtyx.varname = "x", data = mydata) 
```

Note that the method assumes there is only one erroneous variable in the set of independent variables.

### Citation 
```bibtex
@article{satman2015reducing,
  title={Reducing errors-in-variables bias in linear regression using compact genetic algorithms},
  author={Satman, M Hakan and Diyarbakirlioglu, Erkin},
  journal={Journal of Statistical Computation and Simulation},
  volume={85},
  number={16},
  pages={3216--3235},
  year={2015},
  doi={10.1080/00949655.2014.961157}
  publisher={Taylor \& Francis}
}
```