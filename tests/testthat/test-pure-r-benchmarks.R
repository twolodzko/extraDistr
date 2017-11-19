
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
  zb <- sqrt(z/beta)
  bz <- sqrt(beta/z)
  ifelse(x<=mu, 0, 
         (zb+bz)/(2.0*alpha*z) * dnorm((zb-bz)/alpha))
}

# dbvpoisR

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
    1/sigma * (1-xi*z)^(-1-1/xi) * exp(-(1-xi*z)^(-1/xi))
  } else {
    1/sigma * exp(-z) * exp(-exp(-z)) 
  }
}

dgompertzR <- function(x, a, b) {
  ifelse(x<0, 0,
         a*exp(b*x - a/b * (exp(b*x)-1)))
}

dgpdR <- function(x, mu, sigma, xi) {
  z <- (x-mu)/sigma
  if (xi != 0) {
    (1+xi*z)^(-(xi+1)/xi)/sigma
  } else {
    exp(-z)/sigma
  }
}

dgumbelR <- function(x, mu, sigma) {
  z <- (x-mu)/sigma
  1/sigma * exp(-(z+exp(-z)))
}

# dhuberR

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
  ifelse(x<0, 0,
         b*exp(-b*x) * exp(-eta*exp(-b*x)) * exp(1 + eta*(1 - exp(-b*x))))
}

# ====================================================
#                     TESTS
# ====================================================

test_that("Compare PDF's/PMF's to pure-R benchmarks", {
  
  x <- c(-1e5, -100, -10, -5, -1, -0.5, 0.001, 0, 0.001, 0.5, 1, 5, 10, 100, 1e5)
  n <- length(x)
  
  expect_warning(expect_equal(dbbinom(x, 100, 1, 10), dbbinomR(x, 100, 1, 10)))
  expect_warning(expect_equal(dbnbinom(x[1:(n-2)], 100, 1, 10),
                              dbnbinomR(x[1:(n-2)], 100, 1, 10))) # numerical precission in R version
  expect_equal(dbetapr(x, 1, 1, 1), dbetaprR(x, 1, 1, 1))
  # expect_equal(dfatigue(x, 1, 1, 0), dfatigueR(x, 1, 1, 0)) # warning in R version?
  expect_warning(expect_equal(ddlaplace(x, 0, 0.5), ddlaplaceR(x, 0, 0.5)))
  expect_warning(expect_equal(ddweibull(x, 0.5, 1), ddweibullR(x, 0.5, 1)))
  expect_equal(dfrechet(x, 1, 1, 1), dfrechetR(x, 1, 1, 1))
  expect_warning(expect_equal(dgpois(x[1:(n-1)], 1, 1),
                              dgpoisR(x[1:(n-1)], 1, 1))) # numerical precission in R version
  # expect_equal(dgev(x, 1, 1, 1), dgevR(x, 1, 1, 1)) # better rules
  expect_equal(dgompertz(x, 1, 1), dgompertzR(x, 1, 1))
  # expect_equal(dgpd(x, 1, 1, 1), dgpdR(x, 1, 1, 1)) # better rules
  expect_equal(dgumbel(x, 1, 1), dgumbelR(x, 1, 1))
  expect_equal(dinvgamma(x, 1, 1), dinvgammaR(x, 1, 1))
  expect_equal(dlaplace(x, -1, 5), dlaplaceR(x, -1, 5))
  expect_warning(expect_equal(dlgser(x, 0.5), dlgserR(x, 0.5)))
  expect_equal(dlomax(x, 1, 0.5), dlomaxR(x, 1, 0.5))
  expect_equal(dpareto(x, 1, 1), dparetoR(x, 1, 1))
  expect_equal(dpower(x, 1, 1), dpowerR(x, 1, 1))
  expect_equal(drayleigh(x, 1), drayleighR(x, 1))
  # expect_equal(dsgomp(x, 0.5, 1), dsgompR(x, 0.5, 1)) # ???
  
})
