

#' Bivariate Normal distribution
#'
#' Density, distribution function and random generation
#' for the Bivariate Normal distribution.
#'
#' @param x,y                       vectors of quantiles.
#' @param n	                        number of observations. If \code{length(n) > 1},
#'                                  the length is taken to be the number required.
#' @param mu1,mu2,sigma1,sigma2,rho location, scale and correlation parameters.
#' @param log     	                logical; if TRUE, probabilities p are given as log(p).
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
#' For multivariate Normal and t distributions check \pkg{mvtnorm} package.
#'
#' @references
#' Krishnamoorthy, K. (2006). Handbook of Statistical Distributions
#' with Applications. Chapman & Hall/CRC
#'
#' @name BivNormal
#' @aliases BivNormal
#' @aliases dbvnorm
#' @keywords distribution
#'
#' @export

dbvnorm <- function(x, y, mu1 = 0, mu2 = mu1, sigma1 = 1, sigma2 = sigma1, rho = 0, log = FALSE) {
  .Call('extraDistr_cpp_dbnorm', PACKAGE = 'extraDistr', x, y, mu1, mu2, sigma1, sigma2, rho, log)
}


#' @rdname BivNormal
#' @export

rbvnorm <- function(n, mu1 = 0, mu2 = mu1, sigma1 = 1, sigma2 = sigma1, rho = 0) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rbnorm', PACKAGE = 'extraDistr', n, mu1, mu2, sigma1, sigma2, rho)
}

