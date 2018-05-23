

#' Power distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the power distribution.
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
#' hist(x, 100, freq = FALSE)
#' curve(dpower(x, 5, 16), 2, 6, col = "red", add = TRUE, n = 5000)
#' hist(ppower(x, 5, 16))
#' plot(ecdf(x))
#' curve(ppower(x, 5, 16), 2, 6, col = "red", lwd = 2, add = TRUE)
#'
#' @name PowerDist
#' @aliases PowerDist
#' @aliases dpower
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dpower <- function(x, alpha, beta, log = FALSE) {
  cpp_dpower(x, alpha, beta, log[1L])
}


#' @rdname PowerDist
#' @export

ppower <- function(q, alpha, beta, lower.tail = TRUE, log.p = FALSE) {
  cpp_ppower(q, alpha, beta, lower.tail[1L], log.p[1L])
}


#' @rdname PowerDist
#' @export

qpower <- function(p, alpha, beta, lower.tail = TRUE, log.p = FALSE) {
  cpp_qpower(p, alpha, beta, lower.tail[1L], log.p[1L])
}


#' @rdname PowerDist
#' @export

rpower <- function(n, alpha, beta) {
  if (length(n) > 1) n <- length(n)
  cpp_rpower(n, alpha, beta)
}

