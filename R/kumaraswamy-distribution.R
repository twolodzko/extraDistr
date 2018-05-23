

#' Kumaraswamy distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Kumaraswamy distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param a,b             positive valued parameters.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = abx^{a-1} (1-x^a)^{b-1}
#' }{
#' f(x) = a*b*x^(a-1)*(1-x^a)^(b-1)
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = 1-(1-x^a)^b
#' }{
#' F(x) = 1-(1-x^a)^b
#' }
#'
#' Quantile function
#' \deqn{
#' F^{-1}(p) = 1-(1-p^{1/b})^{1/a}
#' }{
#' F^-1(p) = 1-(1-p^(1/b))^(1/a)
#' }
#'
#' @references
#' Jones, M. C. (2009). Kumaraswamy's distribution: A beta-type distribution with
#' some tractability advantages. Statistical Methodology, 6, 70-81.
#'
#' @references
#' Cordeiro, G.M. and de Castro, M. (2009). A new family of generalized distributions.
#' Journal of Statistical Computation & Simulation, 1-17.
#' 
#' @examples 
#' 
#' x <- rkumar(1e5, 5, 16)
#' hist(x, 100, freq = FALSE)
#' curve(dkumar(x, 5, 16), 0, 1, col = "red", add = TRUE)
#' hist(pkumar(x, 5, 16))
#' plot(ecdf(x))
#' curve(pkumar(x, 5, 16), 0, 1, col = "red", lwd = 2, add = TRUE)
#'
#' @name Kumaraswamy
#' @aliases Kumaraswamy
#' @aliases dkumar
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dkumar <- function(x, a = 1, b = 1, log = FALSE) {
  cpp_dkumar(x, a, b, log[1L])
}


#' @rdname Kumaraswamy
#' @export

pkumar <- function(q, a = 1, b = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_pkumar(q, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname Kumaraswamy
#' @export

qkumar <- function(p, a = 1, b = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_qkumar(p, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname Kumaraswamy
#' @export

rkumar <- function(n, a = 1, b = 1) {
  if (length(n) > 1) n <- length(n)
  cpp_rkumar(n, a, b)
}

