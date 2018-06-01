
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

dgumbelR_log <- function(x, mu, sigma) {
  z <- (x-mu)/sigma
  exp(-(z+exp(-z)) - log(sigma))
}

dinvgammaR <- function(x, alpha, beta) {
  ifelse(x<=0, 0, 
         (beta^alpha * x^(-alpha-1) * exp(-beta/x)) / gamma(alpha))
}

dlaplaceR <- function(x, mu, sigma) {
  z <- (x-mu)/sigma
  1/(2*sigma) * exp(-abs(z))
}

dlaplaceR_log <- function(x, mu, sigma) {
  LOG_2F <- 0.6931471805599452862268
  z <- abs(x-mu)/sigma
  exp(-z - LOG_2F - log(sigma))
}

dlgserR <- function(x, theta) {
  ifelse(x<1 | !is_int(x), 0,
         (-1/log(1-theta)*theta^x) / x)
}

dlgserR_log <- function(x, theta) {
  a <- suppressWarnings(-1.0/log1p(-theta))
  ifelse(x<1 | !is_int(x), 0,
         exp(log(a) + (log(theta) * x) - log(x)))
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

