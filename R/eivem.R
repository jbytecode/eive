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
    for (reg.index in 1:yp){
      if (is.null(otherx)) {
        ols <- lm.fit(cbind(1, x_proxy), y[, reg.index])
      }else{
        ols <- lm.fit(cbind(1, x_proxy, otherx), y[, reg.index])
      }
      square_sums <- square_sums + sum(ols$residuals^2)
    }
    return(square_sums)
  }
  g <- cmpfun(f)
  ga <- cga(evalFunc = f, chsize = n * numdummies, popsize = popsize)
  best <- as.integer(ga)
  for (reg.index in 1:yp){
    if (is.null(otherx)) {
      ols_dirty[[reg.index]] <- lm(y[, reg.index] ~ dirtyx)
    }else{
      ols_dirty[[reg.index]] <- lm(y[, reg.index] ~ dirtyx + otherx)
    }
    ols_proxy <- lm(dirtyx ~ matrix(best, nrow = n))
    if (is.null(otherx)) {
      ols_best[[reg.index]] <- lm(y[, reg.index] ~ ols_proxy$fitted.values)
    }else{
      ols_best[[reg.index]] <- lm(y[, reg.index] ~ ols_proxy$fitted.values + otherx)
    }
  } #end of for loop
  result <- list(ols = ols_dirty, eive = ols_best, proxy = ols_proxy)
  return(result)
}
# End of function

