

#' Laplace distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Laplace distribution.
#'
#' @param x,q	            vector of quantiles.
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
#' f(x) = \frac{1}{2\sigma} \exp\left(-\left|\frac{x-\mu}{\sigma}\right|\right)
#' }{
#' f(x) = 1/(2*\sigma) * exp(-|(x-\mu)/\sigma|)
#' }
#'
#' Cumulative distribution function
#' \deqn{ F(x) = \left\{\begin{array}{ll}
#' \frac{1}{2} \exp\left(\frac{x-\mu}{\sigma}\right)     & x < \mu \\
#' 1 - \frac{1}{2} \exp\left(\frac{x-\mu}{\sigma}\right) & x \geq \mu
#' \end{array}\right.
#' }{
#' F(x) = [if x < mu:] 1/2 * exp((x-\mu)/\sigma)
#'        [else:] 1 - 1/2 * exp((x-\mu)/\sigma)
#' }
#'
#' Quantile function
#' \deqn{ F^{-1}(p) = \left\{\begin{array}{ll}
#' \mu + \sigma \log(2p)     & p < 0.5 \\
#' \mu - \sigma \log(2(1-p)) & p \geq 0.5
#' \end{array}\right.
#' }{
#' F^-1(p) = [if p < 0.5:] \mu + \sigma * log(2*p)
#'           [else:] \mu - \sigma * log(2*(1-p))
#' }
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
#' hist(x, 100, freq = FALSE)
#' curve(dlaplace(x, 5, 16), -200, 200, n = 500, col = "red", add = TRUE)
#' hist(plaplace(x, 5, 16))
#' plot(ecdf(x))
#' curve(plaplace(x, 5, 16), -200, 200, n = 500, col = "red", lwd = 2, add = TRUE)
#'
#' @name Laplace
#' @aliases Laplace
#' @aliases dlaplace
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dlaplace <- function(x, mu = 0, sigma = 1, log = FALSE) {
  cpp_dlaplace(x, mu, sigma, log[1L])
}


#' @rdname Laplace
#' @export

plaplace <- function(q, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_plaplace(q, mu, sigma, lower.tail[1L], log.p[1L])
}


#' @rdname Laplace
#' @export

qlaplace <- function(p, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_qlaplace(p, mu, sigma, lower.tail[1L], log.p[1L])
}


#' @rdname Laplace
#' @export

rlaplace <- function(n, mu = 0, sigma = 1) {
  if (length(n) > 1) n <- length(n)
  cpp_rlaplace(n, mu, sigma)
}

