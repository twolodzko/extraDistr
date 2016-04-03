

#' Discrete Weibull distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Discrete Weibull distribution.
#'
#' @param x 	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param q,beta          parameters.
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
#' F^-1(p) = ceiling((\log(1-p)/\log(q))^(1/\beta) - 1)
#' }
#'
#' @references
#' Nakagawa, T. and Osaki, S. (1975). The Discrete Weibull Distribution.
#' IEEE Transactions on Reliability, R-24, 300-301.
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
#' lines(xx, pdweibull(xx, .32, 1), col = "red", lwd = 2)
#'
#' @name DiscrWeibull
#' @aliases DiscrWeibull
#' @aliases ddweibull
#' @keywords distribution
#'
#' @export

ddweibull <- function(x, q, beta, log = FALSE) {
  .Call('extraDistr_cpp_ddweibull', PACKAGE = 'extraDistr', x, q, beta, log)
}


#' @rdname DiscrWeibull
#' @export

pdweibull <- function(x, q, beta, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_pdweibull', PACKAGE = 'extraDistr', x, q, beta, lower.tail, log.p)
}


#' @rdname DiscrWeibull
#' @export

qdweibull <- function(p, q, beta, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qdweibull', PACKAGE = 'extraDistr', p, q, beta, lower.tail, log.p)
}


#' @rdname DiscrWeibull
#' @export

rdweibull <- function(n, q, beta) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rdweibull', PACKAGE = 'extraDistr', n, q, beta)
}

