

#' Generalized Pareto Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Generalized Pareto distribution.
#'
#' @param x 	            vector of quantiles.
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
#' \frac{1}{\sigma} (1+\xi z)^{-(\xi+1)/\xi} & \xi \neq 0 \\
#' \frac{1}{\sigma} \exp(-z)                 & \xi = 0
#' \end{array}\right.
#' }{
#' f(x) = [if \xi != 0:] (1+\xi*z)^{-(\xi+1)/\xi}/sigma
#'        [else:] exp(-z)/sigma
#' }
#'
#' Cumulative distribution function
#' \deqn{ F(x) = \left\{\begin{array}{ll}
#' 1-(1+\xi z)^{-1/\xi} & \xi \neq 0 \\
#' 1-\exp(-z)           & \xi = 0
#' \end{array}\right.
#' }{
#' F(x) = [if \xi != 0:] 1-(1+\xi*z)^{-1/\xi}
#'        [else:] 1-exp(-z)
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
#' where \eqn{ z = \frac{x-\mu}{\sigma} }{ z = (x-\mu)/\sigma }.
#'
#' @references
#' Coles, S. (2001). An Introduction to Statistical Modeling of Extreme Values.
#' Springer.
#' 
#' @examples 
#' 
#' x <- rgpd(1e5, 5, 2, .1)
#' xx <- seq(0, 1000, by = 0.1)
#' hist(x, 100, freq = FALSE, xlim = c(0, 50))
#' lines(xx, dgpd(xx, 5, 2, .1), col = "red")
#' hist(pgpd(x, 5, 2, .1))
#' plot(ecdf(x))
#' lines(xx, pgpd(xx, 5, 2, .1), col = "red", lwd = 2)
#'
#' @name GPD
#' @aliases GPD
#' @aliases dgpd
#' @keywords distribution
#'
#' @export

dgpd <- function(x, mu = 0, sigma = 1, xi = 0, log = FALSE) {
  .Call('extraDistr_cpp_dgpd', PACKAGE = 'extraDistr', x, mu, sigma, xi, log)
}


#' @rdname GPD
#' @export

pgpd <- function(x, mu = 0, sigma = 1, xi = 0, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_pgpd', PACKAGE = 'extraDistr', x, mu, sigma, xi, lower.tail, log.p)
}


#' @rdname GPD
#' @export

qgpd <- function(p, mu = 0, sigma = 1, xi = 0, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qgpd', PACKAGE = 'extraDistr', p, mu, sigma, xi, lower.tail, log.p)
}


#' @rdname GPD
#' @export

rgpd <- function(n, mu = 0, sigma = 1, xi = 0) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rgpd', PACKAGE = 'extraDistr', n, mu, sigma, xi)
}

