

#' Huber Density
#'
#' Density, distribution function, quantile function and random generation
#' for the "Huber density" distribution.
#'
#' @param x 	             vector of quantiles.
#' @param p	               vector of probabilities.
#' @param n	               number of observations. If \code{length(n) > 1},
#'                         the length is taken to be the number required.
#' @param mu,sigma,epsilon location, and scale, and shape parameters.
#' @param log,log.p	       logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	     logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                         otherwise, \eqn{P[X > x]}.
#'
#' @examples 
#' 
#' x <- rhuber(1e5, 5, 2, 3)
#' xx <- seq(-20, 20, by = 0.1)
#' hist(x, 100, freq = FALSE)
#' lines(xx, dhuber(xx, 5, 2, 3), col = "red")
#' hist(phuber(x, 5, 2, 3))
#' plot(ecdf(x))
#' lines(xx, phuber(xx, 5, 2, 3), col = "red", lwd = 2)
#'
#' @name Huber
#' @aliases Huber
#' @aliases dhuber
#' @keywords distribution
#'
#' @export

dhuber <- function(x, mu = 0, sigma = 1, epsilon = 1.345, log = FALSE) {
  .Call('extraDistr_cpp_dhuber', PACKAGE = 'extraDistr', x, mu, sigma, epsilon, log)
}


#' @rdname Huber
#' @export

phuber <- function(x, mu = 0, sigma = 1, epsilon = 1.345, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_phuber', PACKAGE = 'extraDistr', x, mu, sigma, epsilon, lower.tail, log.p)
}


#' @rdname Huber
#' @export

qhuber <- function(p, mu = 0, sigma = 1, epsilon = 1.345, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qhuber', PACKAGE = 'extraDistr', p, mu, sigma, epsilon, lower.tail, log.p)
}


#' @rdname Huber
#' @export

rhuber <- function(n, mu = 0, sigma = 1, epsilon = 1.345) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rhuber', PACKAGE = 'extraDistr', n, mu, sigma, epsilon)
}

