

#' Frechet distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Frechet distribution.
#'
#' @param x 	            vector of quantiles.
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
#' f(x) = \frac{\lambda}{\sigma} z^{-1-\lambda} \exp(-z^{-\lambda})
#' }{
#' f(x) = \lambda/\sigma * z^{-1-\lambda} * exp(-z^-\lambda)
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = \exp(-z^{-\lambda})
#' }{
#' F(x) = exp(-z^-\lambda)
#' }
#'
#' Quantile function
#' \deqn{
#' F^{-1}(p) = \mu + \sigma -\log(p)^{-1/\lambda}
#' }{
#' F^-1(p) = \mu + \sigma * -log(p)^{-1/\lambda}
#' }
#'
#' where \eqn{ z = \frac{x-\mu}{\sigma} }{ z = (x-\mu)/\sigma }.
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
#' @keywords distribution
#'
#' @export

dfrechet <- function(x, lambda = 1, mu = 0, sigma = 1, log = FALSE) {
  .Call('extraDistr_cpp_dfrechet', PACKAGE = 'extraDistr', x, lambda, mu, sigma, log)
}


#' @rdname Frechet
#' @export

pfrechet <- function(x, lambda = 1, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_pfrechet', PACKAGE = 'extraDistr', x, lambda, mu, sigma, lower.tail, log.p)
}


#' @rdname Frechet
#' @export

qfrechet <- function(p, lambda = 1, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qfrechet', PACKAGE = 'extraDistr', p, lambda, mu, sigma, lower.tail, log.p)
}


#' @rdname Frechet
#' @export

rfrechet <- function(n, lambda = 1, mu = 0, sigma = 1) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rfrechet', PACKAGE = 'extraDistr', n, lambda, mu, sigma)
}

