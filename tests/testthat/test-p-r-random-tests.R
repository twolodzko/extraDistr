

# Tests from 'tests/p-r-random-tests.R' from base R

superror <- function(rfoo,pfoo,sample.size,...) {
  x <- rfoo(sample.size,...)
  tx <- table(signif(x, 12)) # such that xi will be sort(unique(x))
  xi <- as.numeric(names(tx))
  f <- pfoo(xi,...)
  fhat <- cumsum(tx)/sample.size
  max(abs(fhat-f))
}

pdkwbound <- function(n,t) 2*exp(-2*n*t*t)

qdkwbound <- function(n,p) sqrt(log(p/2)/(-2*n))

dkwtest <- function(stub = "norm", ...,
                    sample.size = 10000, pthreshold = 0.001,
                    print.result = FALSE, print.detail = FALSE, # don't print by default
                    stop.on.failure = FALSE)                    # don't stop by default
{
  rfoo <- eval(as.name(paste0("r", stub)))
  pfoo <- eval(as.name(paste0("p", stub)))
  s <- superror(rfoo, pfoo, sample.size, ...)
  if (print.result || print.detail) {
    printargs <- substitute(list(...))
    printargs[[1]] <- as.name(stub)
    cat(deparse(printargs))
    if (print.detail)
      cat("\nsupremum error = ",signif(s,2),
          " with p-value=",min(1,round(pdkwbound(sample.size,s),4)),"\n")
  }
  rval <- (s < qdkwbound(sample.size,pthreshold))
  if (print.result)
    cat(c(" FAILED\n"," PASSED\n")[rval+1])
  if (stop.on.failure && !rval)
    stop("dkwtest failed")
  rval
}


test_that("p-r random tests", {
  
  skip_on_cran()
  
  expect_true(dkwtest("bbinom", 1, 1, 1))
  expect_true(dkwtest("bbinom", 10, 1, 1))
  expect_true(dkwtest("bbinom", 10, 100, 1))
  expect_true(dkwtest("bbinom", 10, 1, 100))
  expect_true(dkwtest("bbinom", 100, 1, 1))
  expect_true(dkwtest("bbinom", 100, 100, 1))
  expect_true(dkwtest("bbinom", 100, 1, 100))
  
  expect_true(dkwtest("bern", 0.5))
  
  expect_true(dkwtest("betapr", 1, 1, 1))
  
  expect_true(dkwtest("bhatt", 1, 1, 1))
  
  expect_true(dkwtest("bnbinom", 1, 1, 1))
  expect_true(dkwtest("bnbinom", 10, 1, 1))
  expect_true(dkwtest("bnbinom", 10, 100, 1))
  
  expect_true(dkwtest("cat", c(0.5, 0.5)))
  expect_true(dkwtest("cat", c(1e-6, 1, 2.5, 100)))
  
  expect_true(dkwtest("dlaplace", 1, 0.001))
  expect_true(dkwtest("dlaplace", 0, 0.5))
  expect_true(dkwtest("dlaplace", 0, 0.999))
  
  expect_true(dkwtest("dnorm", 1, 1))
  
  expect_true(dkwtest("dgamma", 9, 1))
  
  expect_true(dkwtest("dunif", 1, 10))
  
  expect_true(dkwtest("dweibull", 0.5, 1))
  expect_true(dkwtest("dweibull", 0.001, 1))
  expect_true(dkwtest("dweibull", 0.999, 1))
  
  expect_true(dkwtest("fatigue", 1, 1))
  
  expect_true(dkwtest("frechet", 1, 1, 1))
  
  expect_true(dkwtest("gev", 1, 1, 1))
  
  expect_true(dkwtest("gompertz", 1, 1))
  
  expect_true(dkwtest("gpd", 1, 1, 1))
  
  expect_true(dkwtest("gpois", 1, 1))
  
  expect_true(dkwtest("gumbel", 1, 1))
  
  expect_true(dkwtest("hcauchy", 1))
  
  expect_true(dkwtest("hnorm", 1))
  
  expect_true(dkwtest("ht", 5, 1))
  
  expect_true(dkwtest("huber", 0, 1, 1))
  
  expect_true(dkwtest("invgamma", 1, 1))
  
  expect_true(dkwtest("invchisq", 1, 1))
  
  expect_true(dkwtest("kumar", 1, 1))
  expect_true(dkwtest("kumar", 100, 1))
  expect_true(dkwtest("kumar", 1, 100))
  expect_true(dkwtest("kumar", 100, 100))
  
  expect_true(dkwtest("laplace", 0, 1))
  expect_true(dkwtest("laplace", 0, 1000))
  
  expect_true(dkwtest("lgser", 0.001))
  expect_true(dkwtest("lgser", 0.5))
  expect_true(dkwtest("lgser", 0.999))
  
  expect_true(dkwtest("lomax", 1, 0.001))
  expect_true(dkwtest("lomax", 1, 0.5))
  expect_true(dkwtest("lomax", 1, 0.999))
  
  expect_true(dkwtest("mixnorm", c(1,2,3), c(1,2,3), c(1/3,1/3,1/3)))
  
  expect_true(dkwtest("mixpois", c(1,2,3), c(1/3,1/3,1/3)))
  
  expect_true(dkwtest("nhyper", 60, 35, 15))
  expect_true(dkwtest("nhyper", 1, 100, 15))
  expect_true(dkwtest("nhyper", 5, 5, 4))
  expect_true(dkwtest("nhyper", 1000, 5, 4))
  
  expect_true(dkwtest("nsbeta", 1, 1, -2, 2))
  
  expect_true(dkwtest("lst", 2, 0, 1))
  
  expect_true(dkwtest("pareto", 1, 1))
  
  expect_true(dkwtest("power", 1, 1))
  expect_true(dkwtest("power", 5, 16))
  
  expect_true(dkwtest("prop", 10, 0.5))
  expect_true(dkwtest("prop", 100, 0.5))
  expect_true(dkwtest("prop", 1000, 0.5))
  expect_true(dkwtest("prop", 10, 0.01))
  expect_true(dkwtest("prop", 10, 0.5))
  expect_true(dkwtest("prop", 10, 0.99))
  
  expect_true(dkwtest("rayleigh", 1))
  
  expect_true(dkwtest("sgomp", 0.4, 1))
  
  expect_true(dkwtest("slash", 1, 1))

  expect_true(dkwtest("tbinom", 200, 0.5, a = 100))
  expect_true(dkwtest("tbinom", 200, 0.5, b = 100))
  
  expect_true(dkwtest("tnorm", 0, 1, -1, 1))
  expect_true(dkwtest("tnorm", 0, 1, -2, 2))
  expect_true(dkwtest("tnorm", 0, 1, 2, Inf))
  expect_true(dkwtest("tnorm", 0, 1, 4, Inf))
  expect_true(dkwtest("tnorm", 0, 1, -Inf, -2))
  expect_true(dkwtest("tnorm", 0, 1, -Inf, -4))
  expect_true(dkwtest("tnorm", 0, 1, -6, -4))
  expect_true(dkwtest("tnorm", 0, 1, 4, 6))
  
  expect_true(dkwtest("tpois", 5, 0))
  expect_true(dkwtest("tpois", 50, 45, 55))
  
  expect_true(dkwtest("triang"))
  expect_true(dkwtest("triang", 0, 1, 0.5))
  
  expect_true(dkwtest("wald", 1, 1))
  
  expect_true(dkwtest("zip", 1, 0.5))
  
  expect_true(dkwtest("zib", 1, 0.5, 0.5))
  
  expect_true(dkwtest("zinb", 1, 0.5, 0.5))
  
})

