

#' Lomax distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Lomax distribution.
#'
#' @param x 	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param lambda,kappa    parameters.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = \frac{\lambda \kappa}{(1+\lambda x)^{\kappa+1}}
#' }{
#' f(x) = \lambda*\kappa / (1+\lambda*x)^(\kappa+1)
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = 1-(1+\lambda x)^{-\kappa}
#' }{
#' F(x) = 1-(1+\lambda*x)^-\kappa
#' }
#'
#' Quantile function
#' \deqn{
#' F^{-1}(p) = \frac{(1-p)^{-1/\kappa} -1}{\lambda}
#' }{
#' F^-1(p) = ((1-p)^(-1/\kappa)-1) / \lambda
#' }
#'
#' @name Lomax
#' @aliases Lomax
#' @aliases dlomax
#' @keywords distribution
#'
#' @export

dlomax <- function(x, lambda, kappa, log = FALSE) {
  .Call('extraDistr_cpp_dlomax', PACKAGE = 'extraDistr', x, lambda, kappa, log)
}


#' @rdname Lomax
#' @export

plomax <- function(x, lambda, kappa, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_plomax', PACKAGE = 'extraDistr', x, lambda, kappa, lower.tail, log.p)
}


#' @rdname Lomax
#' @export

qlomax <- function(p, lambda, kappa, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qlomax', PACKAGE = 'extraDistr', p, lambda, kappa, lower.tail, log.p)
}


#' @rdname Lomax
#' @export

rlomax <- function(n, lambda, kappa) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rlomax', PACKAGE = 'extraDistr', n, lambda, kappa)
}

