

#' Generalized Pareto distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the generalized Pareto distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param mu,sigma,xi	    location, scale, and shape parameters. Scale must be positive.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability density function
#' \deqn{ f(x) = \left\{\begin{array}{ll}
#' \frac{1}{\sigma} \left(1+\xi \frac{x-\mu}{\sigma}\right)^{-(\xi+1)/\xi} & \xi \neq 0 \\
#' \frac{1}{\sigma} \exp\left(-\frac{x-\mu}{\sigma}\right)                 & \xi = 0
#' \end{array}\right.
#' }{
#' f(x) = [if \xi != 0:] (1+\xi*(x-\mu)/\sigma)^{-(\xi+1)/\xi}/\sigma
#'        [else:] exp(-(x-\mu)/\sigma)/\sigma
#' }
#'
#' Cumulative distribution function
#' \deqn{ F(x) = \left\{\begin{array}{ll}
#' 1-\left(1+\xi \frac{x-\mu}{\sigma}\right)^{-1/\xi} & \xi \neq 0 \\
#' 1-\exp\left(-\frac{x-\mu}{\sigma}\right)           & \xi = 0
#' \end{array}\right.
#' }{
#' F(x) = [if \xi != 0:] 1-(1+\xi*(x-\mu)/\sigma)^{-1/\xi}
#'        [else:] 1-exp(-(x-\mu)/\sigma)
#' }
#'
#' Quantile function
#' \deqn{ F^{-1}(x) = \left\{\begin{array}{ll}
#' \mu + \sigma \frac{(1-p)^{-\xi}-1}{\xi} & \xi \neq 0 \\
#' \mu - \sigma \log(1-p)                  & \xi = 0
#' \end{array}\right.
#' }{
#' F^-1(x) = [if \xi != 0:] \mu + \sigma * ((1-p)^{-\xi}-1)/\xi
#'           [else:] \mu - \sigma * log(1-p)
#' }
#'
#' @references
#' Coles, S. (2001). An Introduction to Statistical Modeling of Extreme Values.
#' Springer.
#' 
#' @examples 
#' 
#' x <- rgpd(1e5, 5, 2, .1)
#' hist(x, 100, freq = FALSE, xlim = c(0, 50))
#' curve(dgpd(x, 5, 2, .1), 0, 50, col = "red", add = TRUE, n = 5000)
#' hist(pgpd(x, 5, 2, .1))
#' plot(ecdf(x))
#' curve(pgpd(x, 5, 2, .1), 0, 50, col = "red", lwd = 2, add = TRUE)
#'
#' @name GPD
#' @aliases GPD
#' @aliases dgpd
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dgpd <- function(x, mu = 0, sigma = 1, xi = 0, log = FALSE) {
  cpp_dgpd(x, mu, sigma, xi, log[1L])
}


#' @rdname GPD
#' @export

pgpd <- function(q, mu = 0, sigma = 1, xi = 0, lower.tail = TRUE, log.p = FALSE) {
  cpp_pgpd(q, mu, sigma, xi, lower.tail[1L], log.p[1L])
}


#' @rdname GPD
#' @export

qgpd <- function(p, mu = 0, sigma = 1, xi = 0, lower.tail = TRUE, log.p = FALSE) {
  cpp_qgpd(p, mu, sigma, xi, lower.tail[1L], log.p[1L])
}


#' @rdname GPD
#' @export

rgpd <- function(n, mu = 0, sigma = 1, xi = 0) {
  if (length(n) > 1) n <- length(n)
  cpp_rgpd(n, mu, sigma, xi)
}

