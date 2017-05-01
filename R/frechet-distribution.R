

#' Frechet distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Frechet distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param lambda,sigma,mu shape, scale, and location parameters.
#'                        Scale and shape must be positive.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = \frac{\lambda}{\sigma} \left(\frac{x-\mu}{\sigma}\right)^{-1-\lambda} \exp\left(-\left(\frac{x-\mu}{\sigma}\right)^{-\lambda}\right)
#' }{
#' f(x) = \lambda/\sigma * ((x-\mu)/\sigma)^(-1-\lambda) * exp(-((x-\mu)/\sigma)^-\lambda)
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = \exp\left(-\left(\frac{x-\mu}{\sigma}\right)^{-\lambda}\right)
#' }{
#' F(x) = exp(-((x-\mu)/\sigma)^-\lambda)
#' }
#'
#' Quantile function
#' \deqn{
#' F^{-1}(p) = \mu + \sigma -\log(p)^{-1/\lambda}
#' }{
#' F^-1(p) = \mu + \sigma * -log(p)^{-1/\lambda}
#' }
#'
#' @references
#' Bury, K. (1999). Statistical Distributions in Engineering.
#' Cambridge University Press.
#' 
#' @examples 
#' 
#' x <- rfrechet(1e5, 5, 2, 1.5)
#' xx <- seq(0, 1000, by = 0.1)
#' hist(x, 200, freq = FALSE)
#' lines(xx, dfrechet(xx, 5, 2, 1.5), col = "red") 
#' hist(pfrechet(x, 5, 2, 1.5))
#' plot(ecdf(x))
#' lines(xx, pfrechet(xx, 5, 2, 1.5), col = "red", lwd = 2)
#'
#' @name Frechet
#' @aliases Frechet
#' @aliases dfrechet
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dfrechet <- function(x, lambda = 1, mu = 0, sigma = 1, log = FALSE) {
  cpp_dfrechet(x, lambda, mu, sigma, log[1L])
}


#' @rdname Frechet
#' @export

pfrechet <- function(q, lambda = 1, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_pfrechet(q, lambda, mu, sigma, lower.tail[1L], log.p[1L])
}


#' @rdname Frechet
#' @export

qfrechet <- function(p, lambda = 1, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_qfrechet(p, lambda, mu, sigma, lower.tail[1L], log.p[1L])
}


#' @rdname Frechet
#' @export

rfrechet <- function(n, lambda = 1, mu = 0, sigma = 1) {
  if (length(n) > 1) n <- length(n)
  cpp_rfrechet(n, lambda, mu, sigma)
}

