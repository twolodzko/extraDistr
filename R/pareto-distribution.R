

#' Pareto distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Pareto distribution.
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
#' f(x) = \frac{ab^a}{x^{a+1}}
#' }{
#' f(x) = (a*b^a) / x^(a+1)
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = 1 - \left(\frac{b}{x}\right)^a
#' }{
#' F(x) = 1 - (b/x)^a
#' }
#'
#' Quantile function
#' \deqn{
#' F^{-1}(p) = \frac{b}{(1-p)^{1-a}}
#' }{
#' F^-1(p) = b/(1-p)^(1-a)
#' }
#'
#' @references
#' Krishnamoorthy, K. (2006). Handbook of Statistical Distributions
#' with Applications. Chapman & Hall/CRC
#' 
#' @examples 
#' 
#' x <- rpareto(1e5, 5, 16)
#' hist(x, 100, freq = FALSE)
#' curve(dpareto(x, 5, 16), 0, 200, col = "red", add = TRUE)
#' hist(ppareto(x, 5, 16))
#' plot(ecdf(x))
#' curve(ppareto(x, 5, 16), 0, 200, col = "red", lwd = 2, add = TRUE)
#'
#' @name Pareto
#' @aliases Pareto
#' @aliases dpareto
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dpareto <- function(x, a = 1, b = 1, log = FALSE) {
  cpp_dpareto(x, a, b, log[1L])
}


#' @rdname Pareto
#' @export

ppareto <- function(q, a = 1, b = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_ppareto(q, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname Pareto
#' @export

qpareto <- function(p, a = 1, b = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_qpareto(p, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname Pareto
#' @export

rpareto <- function(n, a = 1, b = 1) {
  if (length(n) > 1) n <- length(n)
  cpp_rpareto(n, a, b)
}

