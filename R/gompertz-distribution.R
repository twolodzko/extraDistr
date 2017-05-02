

#' Gompertz distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Gompertz distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param a,b             positive valued scale and location parameters.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = a \exp\left(bx - \frac{a}{b} (\exp(bx)-1)\right)
#' }{
#' f(x) = a*exp(b*x - a/b * (exp(b*x)-1))
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = 1-\exp\left(-\frac{a}{b} (\exp(bx)-1)\right)
#' }{
#' F(x) = 1-exp(-a/b * (exp(b*x)-1))
#' }
#'
#' Quantile function
#' \deqn{
#' F^{-1}(p) = \frac{1}{b} \log\left(1-\frac{b}{a}\log(1-p)\right)
#' }{
#' F^-1(p) = 1/b * log(1 - b/a * log(1-p))
#' }
#'
#' @references
#' Lenart, A. (2012). The Gompertz distribution and Maximum Likelihood Estimation
#' of its parameters - a revision. MPIDR WORKING PAPER WP 2012-008.
#' \url{http://www.demogr.mpg.de/papers/working/wp-2012-008.pdf}
#' 
#' @examples 
#' 
#' x <- rgompertz(1e5, 5, 2)
#' hist(x, 100, freq = FALSE)
#' curve(dgompertz(x, 5, 2), 0, 1, col = "red", add = TRUE)
#' hist(pgompertz(x, 5, 2))
#' plot(ecdf(x))
#' curve(pgompertz(x, 5, 2), 0, 1, col = "red", lwd = 2, add = TRUE)
#'
#' @name Gompertz
#' @aliases Gompertz
#' @aliases dgompertz
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dgompertz <- function(x, a = 1, b = 1, log = FALSE) {
  cpp_dgompertz(x, a, b, log[1L])
}


#' @rdname Gompertz
#' @export

pgompertz <- function(q, a = 1, b = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_pgompertz(q, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname Gompertz
#' @export

qgompertz <- function(p, a = 1, b = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_qgompertz(p, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname Gompertz
#' @export

rgompertz <- function(n, a = 1, b = 1) {
  if (length(n) > 1) n <- length(n)
  cpp_rgompertz(n, a, b)
}

