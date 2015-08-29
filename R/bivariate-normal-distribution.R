

#' Bivariate Normal distribution
#'
#' Density, distribution function and random generation
#' for the Bivariate Normal distribution.
#'
#' @param x                         matrix of quantiles.
#' @param n	                        number of observations. If \code{length(n) > 1},
#'                                  the length is taken to be the number required.
#' @param mu1,mu2,sigma1,sigma2,rho location, scale and correlation parameters.
#' @param log,log.p	                logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	              logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                                  otherwise, \eqn{P[X > x]}.
#' @param nsim                      number of samples in Monte Carlo simulation for calculating
#'                                  cumulative distribution function; the higher is more precise
#'                                  but slower.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = \frac{1}{2\pi\sqrt{1-\rho^2}\sigma_1\sigma_2}
#'        \exp\left(-\frac{1}{2(1-\rho^2)} (z_1^2 - 2\rho z_1 z_2 + z_2^2)\right)
#' }{
#' f(x) = 1/(2*pi*sqrt(1-rho^2)*sigma1*sigma2) * exp(-(1/(2*(1-rho^2)*(z1^2 - 2*rho*z1*z2 * z2^2))))
#' }
#'
#' where \eqn{
#' z_1 = \frac{x_1 - \mu_1}{\sigma_1}
#' }{
#' z1 = (x1 - \mu1)/\sigma1
#' } and \eqn{
#' z_2 = \frac{x_2 - \mu_2}{\sigma_2}
#' }{
#' z2 = (x2 - \mu2)/\sigma2
#' }.
#'
#' Cumulative distribution function is approximated using Monte Carlo simulation.
#' For multivariate Normal and t distributions check \pkg{mvtnorm} package.
#'
#' @references
#' Krishnamoorthy, K. (2006). Handbook of Statistical Distributions
#' with Applications. Chapman & Hall/CRC
#'
#' @name BivNormal
#' @aliases BivNormal
#' @aliases dbnorm
#' @keywords distribution
#'
#' @export

dbnorm <- function(x, mu1 = 0, mu2 = mu1, sigma1 = 1, sigma2 = sigma1, rho = 0, log = FALSE) {
  .Call('extraDistr_cpp_dbnorm', PACKAGE = 'extraDistr', x, mu1, mu2, sigma1, sigma2, rho, log)
}


#' @rdname BivNormal
#' @export

pbnorm <- function(x, mu1 = 0, mu2 = mu1, sigma1 = 1, sigma2 = sigma1, rho = 0,
                   lower.tail = TRUE, log.p = FALSE, nsim = 10000L) {
  .Call('extraDistr_cpp_pbnorm', PACKAGE = 'extraDistr', x, mu1, mu2, sigma1, sigma2, rho, lower.tail, log.p, nsim)
}


#' @rdname BivNormal
#' @export

rbnorm <- function(n, mu1 = 0, mu2 = mu1, sigma1 = 1, sigma2 = sigma1, rho = 0) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rbnorm', PACKAGE = 'extraDistr', n, mu1, mu2, sigma1, sigma2, rho)
}

