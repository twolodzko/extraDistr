

#' Discrete Weibull distribution (type I)
#'
#' Density, distribution function, quantile function and random generation
#' for the discrete Weibull (type I) distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param shape1,shape2   parameters (named q, \eqn{\beta}). Values of \code{shape2}
#'                        need to be positive and \code{0 < shape1 < 1}.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability mass function
#' \deqn{
#' f(x) = q^{x^\beta} - q^{(x+1)^\beta}
#' }{
#' f(x) = q^x^\beta - q^(x+1)^\beta
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = 1-q^{(x+1)^\beta}
#' }{
#' F(x) = 1-q^(x+1)^\beta
#' }
#'
#' Quantile function
#' \deqn{
#' F^{-1}(p) = \left \lceil{\left(\frac{\log(1-p)}{\log(q)}\right)^{1/\beta} - 1}\right \rceil
#' }{
#' F^-1(p) = ceiling((log(1-p)/log(q))^(1/\beta) - 1)
#' }
#'
#' @references
#' Nakagawa, T. and Osaki, S. (1975). The Discrete Weibull Distribution.
#' IEEE Transactions on Reliability, R-24, 300-301.
#' 
#' @references 
#' Kulasekera, K.B. (1994).
#' Approximate MLE's of the parameters of a discrete Weibull distribution
#' with type I censored data.
#' Microelectronics Reliability, 34(7), 1185-1188.
#' 
#' @references 
#' Khan, M.A., Khalique, A. and Abouammoh, A.M. (1989).
#' On estimating parameters in a discrete Weibull distribution.
#' IEEE Transactions on Reliability, 38(3), 348-350.
#' 
#' @seealso \code{\link[stats]{Weibull}}
#' 
#' @examples 
#' 
#' x <- rdweibull(1e5, 0.32, 1)
#' xx <- seq(-2, 100, by = 1)
#' plot(prop.table(table(x)), type = "h")
#' lines(xx, ddweibull(xx, .32, 1), col = "red")
#' 
#' # Notice: distribution of F(X) is far from uniform:
#' hist(pdweibull(x, .32, 1), 50)
#' 
#' plot(ecdf(x))
#' lines(xx, pdweibull(xx, .32, 1), col = "red", lwd = 2, type = "s")
#'
#' @name DiscreteWeibull
#' @aliases DiscreteWeibull
#' @aliases ddweibull
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Discrete
#'
#' @export

ddweibull <- function(x, shape1, shape2, log = FALSE) {
  cpp_ddweibull(x, shape1, shape2, log[1L])
}


#' @rdname DiscreteWeibull
#' @export

pdweibull <- function(q, shape1, shape2, lower.tail = TRUE, log.p = FALSE) {
  cpp_pdweibull(q, shape1, shape2, lower.tail[1L], log.p[1L])
}


#' @rdname DiscreteWeibull
#' @export

qdweibull <- function(p, shape1, shape2, lower.tail = TRUE, log.p = FALSE) {
  cpp_qdweibull(p, shape1, shape2, lower.tail[1L], log.p[1L])
}


#' @rdname DiscreteWeibull
#' @export

rdweibull <- function(n, shape1, shape2) {
  if (length(n) > 1) n <- length(n)
  cpp_rdweibull(n, shape1, shape2)
}

