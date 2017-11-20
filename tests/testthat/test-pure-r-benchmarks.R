
is_int <- function(x) x %% 1 == 0

dbbinomR <- function(k, n, alpha, beta) {
  ifelse(k<0 | k>n | !is_int(k), 0,
         choose(n, k) * (beta(k+alpha, n-k+beta)) / (beta(alpha, beta)))
}

dbnbinomR <- function(k, r, alpha, beta) {
  ifelse(k<0 | !is_int(k), 0,
         gamma(r+k)/(factorial(k) * gamma(r)) * beta(alpha+r, beta+k)/beta(alpha, beta))
}

dbetaprR <- function(x, alpha, beta, sigma) {
  z <- x/sigma
  ifelse(x<=0, 0,
         z^(alpha-1.0) * (z+1.0)^(-alpha-beta) / beta(alpha, beta) / sigma)
}

dfatigueR <- function(x, alpha, beta, mu) {
  z <- x-mu
  zb <- suppressWarnings(sqrt(z/beta)) # warns when x <= mu
  bz <- suppressWarnings(sqrt(beta/z)) # warns when x <= mu
  ifelse(x<=mu, 0, 
         (zb+bz)/(2.0*alpha*z) * dnorm((zb-bz)/alpha))
}

ddlaplaceR <- function(x, mu, p) {
  ifelse(!is_int(x), 0,
         (1.0-p)/(1.0+p) * p^abs(x-mu))
}

ddweibullR <- function(x, q, beta) {
  ifelse(x<0 | !is_int(x), 0, 
         q^x^beta - q^(x+1)^beta)
}

dfrechetR <- function(x, lambda, mu, sigma) {
  z <- (x-mu)/sigma
  ifelse(x<=mu, 0,
         lambda/sigma * z^(-1-lambda) * exp(-z^(-lambda)))
}

dgpoisR <- function(x, alpha, beta) {
  ifelse(x<0 | !is_int(x), 0,
         (gamma(x+beta)*alpha^x)/(gamma(beta)*(1+alpha)^(beta+x)*factorial(x)))
}

dgevR <- function(x, mu, sigma, xi) {
  z <- (x-mu)/sigma
  if (xi != 0) {
    ifelse((1.0+xi*z) <= 0.0, 0, 
           1/sigma * (1+xi*z)^(-1/xi-1) * exp(-(1+xi*z)^(-1/xi)))
  } else {
    ifelse((1.0+xi*z) <= 0.0, 0, 
           1/sigma * exp(-z) * exp(-exp(-z)))
  }
}

dgompertzR <- function(x, a, b) {
  ifelse(x<0, 0,
         a*exp(b*x - a/b * (exp(b*x)-1)))
}

dgpdR <- function(x, mu, sigma, xi) {
  z <- (x-mu)/sigma
  if (xi != 0) {
    ifelse((1.0+xi*z) <= 0.0 | z <= 0, 0,
           (1+xi*z)^(-(xi+1)/xi)/sigma)
  } else {
    ifelse((1.0+xi*z) <= 0.0 | z <= 0, 0,
           exp(-z)/sigma)
  }
}

dgumbelR <- function(x, mu, sigma) {
  z <- (x-mu)/sigma
  1/sigma * exp(-(z+exp(-z)))
}

dinvgammaR <- function(x, alpha, beta) {
  ifelse(x<=0, 0, 
         (x^(-alpha-1) * exp(-1/(beta*x))) / (gamma(alpha)*beta^alpha))
}

dlaplaceR <- function(x, mu, sigma) {
  z <- (x-mu)/sigma
  1/(2*sigma) * exp(-abs(z))
}

dlgserR <- function(x, theta) {
  ifelse(x<1 | !is_int(x), 0, 
         (-1/log(1-theta)*theta^x) / x)
}

dlomaxR <- function(x, lambda, kappa) {
  ifelse(x<=0, 0, 
         lambda*kappa / (1+lambda*x)^(kappa+1))
}

dparetoR <- function(x, a, b) {
  ifelse(x<b, 0,
         (a*b^a) / x^(a+1))
}

dpowerR <- function(x, alpha, beta) {
  ifelse(x<=0 | x>=alpha, 0,
         (beta*x^(beta-1)) / (alpha^beta))
}

drayleighR <- function(x, sigma) {
  ifelse(x<=0, 0,
         x/sigma^2 * exp(-(x^2 / 2*sigma^2)))
}

dsgompR <- function(x, b, eta) {
  ebx <- exp(-b*x)
  ifelse(x<0, 0,
         b*ebx * exp(-eta*ebx) * (1 + eta*(1 - ebx)))
}

# ====================================================
#                     TESTS
# ====================================================

x <- c(-1e5, -100, -10, -5.5, -5, -1.01, -1, -0.5, 0.001, 0, 0.001, 0.5, 1, 1.01, 5, 5.5, 10, 100, 1e5)
n <- length(x)


test_that("Compare PDF's/PMF's to pure-R benchmarks", {
  
  expect_warning(expect_equal(dbbinom(x, 100, 1, 10), dbbinomR(x, 100, 1, 10)))
  expect_warning(expect_equal(dbnbinom(x[1:(n-2)], 100, 1, 10),
                              dbnbinomR(x[1:(n-2)], 100, 1, 10))) # numerical precission in R version
  expect_equal(dbetapr(x, 1, 1, 1), dbetaprR(x, 1, 1, 1))
  expect_equal(dfatigue(x, 1, 1, 0), dfatigueR(x, 1, 1, 0))
  expect_warning(expect_equal(ddlaplace(x, 0, 0.5), ddlaplaceR(x, 0, 0.5)))
  expect_warning(expect_equal(ddweibull(x, 0.5, 1), ddweibullR(x, 0.5, 1)))
  expect_equal(dfrechet(x, 1, 1, 1), dfrechetR(x, 1, 1, 1))
  expect_warning(expect_equal(dgpois(x[1:(n-1)], 1, 1),
                              dgpoisR(x[1:(n-1)], 1, 1))) # numerical precission in R version
  expect_equal(dgev(x, 1, 1, 1), dgevR(x, 1, 1, 1))
  expect_equal(dgompertz(x, 1, 1), dgompertzR(x, 1, 1))
  expect_equal(dgpd(x, 1, 1, 1), dgpdR(x, 1, 1, 1))
  expect_equal(dgumbel(x, 1, 1), dgumbelR(x, 1, 1))
  expect_equal(dinvgamma(x, 1, 1), dinvgammaR(x, 1, 1))
  expect_equal(dlaplace(x, -1, 5), dlaplaceR(x, -1, 5))
  expect_warning(expect_equal(dlgser(x, 0.5), dlgserR(x, 0.5)))
  expect_equal(dlomax(x, 1, 0.5), dlomaxR(x, 1, 0.5))
  expect_equal(dpareto(x, 1, 1), dparetoR(x, 1, 1))
  expect_equal(dpower(x, 1, 1), dpowerR(x, 1, 1))
  expect_equal(drayleigh(x, 1), drayleighR(x, 1))
  expect_equal(dsgomp(x, 0.5, 1), dsgompR(x, 0.5, 1))
  
})

test_that("Compare dhuber to hoa implementation", {
  
  skip_if_not_installed("hoa")

  expect_equal(dhuber(x), hoa::dHuber(x))
  expect_equal(phuber(x), hoa::pHuber(x))
  
})

test_that("Compare GEV and GPD to evd implementation", {
  
  skip_if_not_installed("evd")

  expect_equal(dgev(x), evd::dgev(x))
  expect_equal(pgev(x), evd::pgev(x))
  
  expect_equal(dgpd(x), evd::dgpd(x))
  expect_equal(pgpd(x), evd::pgpd(x))
  
})


# ====================================================
#  Compare dbvpois with code by Karlis and Ntzoufras
# ====================================================

# Below code was taken from bivpois package ver 0.50-3 (CRAN archive)

pbivpois <- function(x, y=NULL, a, b, c, log=FALSE) {
  # ------------------------------------------------------------------------------
  # Karlis and Ntzoufras (2003, 2004)
  # EM algorithms for Bivariate Poisson Models
  # ------------------------------------------------------------------------------
  # x      : matrix or vector of length n
  # y      : vector of length n. If x is matrix then it is not used
  # lambda : parameters of the bivariate poisson distribution
  # log    : argument controlling the calculation of the log-probability or the 
  #          probability function. 
  # ------------------------------------------------------------------------------
  #	
  
  lambda <- c(a, b, c) # edited here
  
  if ( is.matrix(x) ) {
    var1<-x[,1]
    var2<-x[,2]
  }
  else if (is.vector(x)&is.vector(y)){
    if (length(x)==length(y)){
      var1<-x
      var2<-y
    }
    else{
      stop('lengths of x and y are not equal')
    }	
  }
  else{
    stop('x is not a matrix or x and y are not vectors')
  }
  n <- length(var1)
  logbp<-vector(length=n)
  #
  for (k in 1:n){
    x0<-var1[k]
    y0<-var2[k]
    xymin<-min( x0,y0 )
    lambdaratio<-lambda[3]/(lambda[1]*lambda[2])
    #	
    i<-0:xymin
    sums<- -lgamma(var1[k]-i+1)-lgamma(i+1)-lgamma(var2[k]-i+1)+i*log(lambdaratio)
    maxsums <- max(sums)
    sums<- sums - maxsums
    logsummation<- log( sum(exp(sums)) ) + maxsums 
    logbp[k]<- -sum(lambda) + var1[k] * log( lambda[1] ) + var2[k] * log( lambda[2] ) + logsummation 
  }
  if (log) { result<-    logbp }
  else     { result<-exp(logbp)  }
  result
  #	end of function bivpois
}


test_that("Bivariate Poisson distribution agrees with code by Karlis and Ntzoufras", {
  
  x <- rbvpois(1000, 7, 8, 5)

  expect_equal(
    dbvpois(x[,1], x[,2], 7, 8, 5),
    pbivpois(x[,1], x[,2], 7, 8, 5)
  )
  
})
