

#' Generalized extreme value distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the generalized extreme value distribution.
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
#' \frac{1}{\sigma} (1-\xi z)^{-1-1/\xi} \exp(-(1-\xi z)^{-1/\xi}) & \xi \neq 0 \\
#' \frac{1}{\sigma} \exp(-z) \exp(-\exp(-z))                       & \xi = 0
#' \end{array}\right.
#' }{
#' f(x) = [if \xi != 0:] 1/\sigma * (1-\xi*z)^{-1-1/\xi} * exp(-(1-\xi*z)^{-1/\xi})
#' [else:] 1/\sigma * exp(-z) * exp(-exp(-z))
#' }
#'
#' Cumulative distribution function
#' \deqn{ F(x) = \left\{\begin{array}{ll}
#' \exp(-(1+\xi z)^{1/\xi}) & \xi \neq 0 \\
#' \exp(-\exp(-z))          & \xi = 0
#' \end{array}\right.
#' }{
#' F(x) = [if \xi != 0:] exp(-(1+\xi*z)^{1/\xi})
#' [else:] exp(-exp(-z))
#' }
#'
#' Quantile function
#' \deqn{ F^{-1}(p) = \left\{\begin{array}{ll}
#' \mu - \frac{\sigma}{\xi} (1 - (-\log(p))^\xi)  & \xi \neq 0 \\
#' \mu - \sigma \log(-\log(p))                    & \xi = 0
#' \end{array}\right.
#' }{
#' F^-1(p) = [if \xi != 0:] \mu - \sigma/\xi * (1 - (-log(p))^\xi)
#'           [else:] \mu - \sigma * log(-log(p))
#' }
#'
#' where \eqn{ z = \frac{x-\mu}{\sigma} }{ z = (x-\mu)/\sigma }.
#'
#' @references
#' Coles, S. (2001). An Introduction to Statistical Modeling of Extreme Values.
#' Springer.
#'
#' @examples 
#' 
#' x <- rgev(1e5, 5, 2, .5)
#' xx <- seq(0, 1000, by = 0.1)
#' hist(x, 1000, freq = FALSE, xlim = c(0, 50))
#' lines(xx, dgev(xx, 5, 2, .5), col = "red")
#' hist(pgev(x, 5, 2, .5))
#' plot(ecdf(x))
#' lines(xx, pgev(xx, 5, 2, .5), col = "red", lwd = 2)
#'
#' @name GEV
#' @aliases GEV
#' @aliases dgev
#' @keywords distribution
#'
#' @export

dgev <- function(x, mu = 0, sigma = 1, xi = 0, log = FALSE) {
  .Call('extraDistr_cpp_dgev', PACKAGE = 'extraDistr', x, mu, sigma, xi, log)
}


#' @rdname GEV
#' @export

pgev <- function(q, mu = 0, sigma = 1, xi = 0, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_pgev', PACKAGE = 'extraDistr', q, mu, sigma, xi, lower.tail, log.p)
}


#' @rdname GEV
#' @export

qgev <- function(p, mu = 0, sigma = 1, xi = 0, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qgev', PACKAGE = 'extraDistr', p, mu, sigma, xi, lower.tail, log.p)
}


#' @rdname GEV
#' @export

rgev <- function(n, mu = 0, sigma = 1, xi = 0) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rgev', PACKAGE = 'extraDistr', n, mu, sigma, xi)
}

