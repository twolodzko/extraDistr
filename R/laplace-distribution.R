

#' Laplace distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Laplace distribution.
#'
#' @param x 	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param mu,sigma	      location and scale parameters. Scale must be positive.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = \frac{1}{2\sigma} \exp(-|z|)
#' }{
#' f(x) = 1/(2*\sigma) * exp(-|z|)
#' }
#'
#' Cumulative distribution function
#' \deqn{ F(x) = \left\{\begin{array}{ll}
#' \frac{1}{2} \exp(z)     & x < \mu \\
#' 1 - \frac{1}{2} \exp(z) & x \geq \mu
#' \end{array}\right.
#' }{
#' F(x) = [if x < mu:] 1/2 * exp(z)
#'        [else:] 1 - 1/2 * exp(z)
#' }
#'
#' Quantile function
#' \deqn{ F^{-1}(p) = \left\{\begin{array}{ll}
#' \mu + \sigma \log(2p)     & p < 0.5 \\
#' \mu + \sigma \log(2(1-p)) & p \geq 0.5
#' \end{array}\right.
#' }{
#' F^-1(p) = [if p < 0.5:] \mu + \sigma * log(2*p)
#'           [else:] \mu + \sigma * log(2*(1-p))
#' }
#'
#' where \eqn{ z = \frac{x-\mu}{\sigma} }{ z = (x-\mu)/\sigma }.
#'
#' @references
#' Krishnamoorthy, K. (2006). Handbook of Statistical Distributions
#' with Applications. Chapman & Hall/CRC
#'
#' @references
#' Forbes, C., Evans, M. Hastings, N., & Peacock, B. (2011).
#' Statistical Distributions. John Wiley & Sons.
#' 
#' @examples 
#' 
#' x <- rlaplace(1e5, 5, 16)
#' xx <- seq(-100, 100, by = 0.01)
#' hist(x, 100, freq = FALSE)
#' lines(xx, dlaplace(xx, 5, 16), col = "red")
#' hist(plaplace(x, 5, 16))
#' plot(ecdf(x))
#' lines(xx, plaplace(xx, 5, 16), col = "red", lwd = 2)
#'
#' @name Laplace
#' @aliases Laplace
#' @aliases dlaplace
#' @keywords distribution
#'
#' @export

dlaplace <- function(x, mu = 0, sigma = 1, log = FALSE) {
  .Call('extraDistr_cpp_dlaplace', PACKAGE = 'extraDistr', x, mu, sigma, log)
}


#' @rdname Laplace
#' @export

plaplace <- function(x, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_plaplace', PACKAGE = 'extraDistr', x, mu, sigma, lower.tail, log.p)
}


#' @rdname Laplace
#' @export

qlaplace <- function(p, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qlaplace', PACKAGE = 'extraDistr', p, mu, sigma, lower.tail, log.p)
}


#' @rdname Laplace
#' @export

rlaplace <- function(n, mu = 0, sigma = 1) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rlaplace', PACKAGE = 'extraDistr', n, mu, sigma)
}

