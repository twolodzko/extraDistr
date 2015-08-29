

#' Beta-Binomial distribution
#'
#' Density, distribution function and random generation
#' for the Beta-binomial distribution.
#'
#' @param x 	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param size,alpha,beta parameters.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability mass function
#' \deqn{
#' f(x) = {n \choose x} \frac{\mathrm{B}(x+\alpha, n-x+\beta)}{\mathrm{B}(\alpha, \beta)}
#' }{
#' f(x) = choose(n, x) * (beta(x+alpha, n-x+beta)) / (beta(alpha, beta))
#' }
#'
#' @seealso \code{\link{Beta}}, \code{\link{Binomial}}
#'
#' @name BetaBinom
#' @aliases BetaBinom
#' @aliases dbbinom
#' @keywords distribution
#'
#' @export

dbbinom <- function(x, size, alpha = 1, beta = 1, log = FALSE) {
  .Call('extraDistr_cpp_dbbinom', PACKAGE = 'extraDistr', x, size, alpha, beta, log)
}


#' @rdname BetaBinom
#' @export

pbbinom <- function(x, size, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_pbbinom', PACKAGE = 'extraDistr', x, size, alpha, beta, lower.tail, log.p)
}


#' @rdname BetaBinom
#' @export

rbbinom <- function(n, size, alpha = 1, beta = 1) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rbbinom', PACKAGE = 'extraDistr', n, size, alpha, beta)
}

