

#' Half-t distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the half-t distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param nu,sigma	      positive valued degrees of freedom and scale parameters.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'                        
#' @details
#' If \eqn{X} follows t distribution parametrized by degrees of freedom \eqn{\nu}
#' and scale \eqn{\sigma}, then \eqn{|X|} follows half-t distribution parametrized
#' by degrees of freedom \eqn{\nu} and scale \eqn{\sigma}.
#'                        
#' @references
#' Gelman, A. (2006). Prior distributions for variance parameters in hierarchical
#' models (comment on article by Browne and Draper).
#' Bayesian analysis, 1(3), 515-534.
#' 
#' @references 
#' Jacob, E. and Jayakumar, K. (2012).
#' On Half-Cauchy Distribution and Process.
#' International Journal of Statistika and Mathematika, 3(2), 77-81.
#' 
#' @seealso \code{\link{HalfNormal}}, \code{\link{HalfCauchy}}
#' 
#' @examples 
#' 
#' x <- rht(1e5, 2, 2)
#' xx <- seq(-1, 100, by = 0.01)
#' hist(x, 500, freq = FALSE, xlim = c(0, 100))
#' lines(xx, dht(xx, 2, 2), col = "red")
#' hist(pht(x, 2, 2))
#' plot(ecdf(x), xlim = c(0, 100))
#' lines(xx, pht(xx, 2, 2), col = "red", lwd = 2)
#'
#' @name HalfT
#' @aliases HalfT
#' @aliases dht
#' @keywords distribution
#'
#' @export

dht <- function(x, nu, sigma = 1, log = FALSE) {
  .Call('extraDistr_cpp_dhalf', PACKAGE = 'extraDistr', x, nu, sigma, log)
}


#' @rdname HalfT
#' @export

pht <- function(q, nu, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_phalf', PACKAGE = 'extraDistr', q, nu, sigma, lower.tail, log.p)
}


#' @rdname HalfT
#' @export

qht <- function(p, nu, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qhalf', PACKAGE = 'extraDistr', p, nu, sigma, lower.tail, log.p)
}


#' @rdname HalfT
#' @export

rht <- function(n, nu, sigma = 1) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rhalf', PACKAGE = 'extraDistr', n, nu, sigma)
}

