

#' Power distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Pareto distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param alpha,beta      parameters.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = \frac{\beta x^{\beta-1}}{\alpha^\beta}
#' }{
#' f(x) = (\beta*x^(\beta-1)) / (\alpha^\beta)
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = \frac{x^\beta}{\alpha^\beta}
#' }{
#' F(x) = x^\beta / \alpha^\beta
#' }
#'
#' Quantile function
#' \deqn{
#' F^{-1}(p) = \alpha p^{1/\beta}
#' }{
#' F^-1(p) = \alpha * p^(1/\beta)
#' }
#' 
#' @examples 
#' 
#' x <- rpower(1e5, 5, 16)
#' xx <- seq(-100, 100, by = 0.001)
#' hist(x, 100, freq = FALSE)
#' lines(xx, dpower(xx, 5, 16), col = "red")
#' hist(ppower(x, 5, 16))
#' plot(ecdf(x))
#' lines(xx, ppower(xx, 5, 16), col = "red", lwd = 2)
#'
#' @name PowerDist
#' @aliases PowerDist
#' @aliases dpower
#' @keywords distribution
#'
#' @export

dpower <- function(x, alpha, beta, log = FALSE) {
  .Call('extraDistr_cpp_dpower', PACKAGE = 'extraDistr', x, alpha, beta, log)
}


#' @rdname PowerDist
#' @export

ppower <- function(q, alpha, beta, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_ppower', PACKAGE = 'extraDistr', q, alpha, beta, lower.tail, log.p)
}


#' @rdname PowerDist
#' @export

qpower <- function(p, alpha, beta, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qpower', PACKAGE = 'extraDistr', p, alpha, beta, lower.tail, log.p)
}


#' @rdname PowerDist
#' @export

rpower <- function(n, alpha, beta) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rpower', PACKAGE = 'extraDistr', n, alpha, beta)
}

