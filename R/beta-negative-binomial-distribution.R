

#' Beta-Negative Binomial distribution
#'
#' Probability mass function and random generation
#' for the Beta-binomial distribution.
#'
#' @param x 	            vector of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param size,alpha,beta parameters.
#' @param log     	      logical; if TRUE, probabilities p are given as log(p).
#'
#' @details
#'
#' Probability mass function
#' \deqn{
#' f(x) = \frac{\Gamma(r+k)}{k! \Gamma(r)}
#'        \frac{\mathrm{B}(\alpha+r, \beta+k)}{\mathrm{B}(\alpha, \beta)}
#' }{
#' f(x) = gamma(r+k)/(k! gamma(r)) * beta(alpha+r, beta+k)/beta(alpha, beta)
#' }
#'
#' @seealso \code{\link{Beta}}, \code{\link{NegBinomial}}
#'
#' @name BetaNegBinom
#' @aliases BetaNegBinom
#' @aliases dbnbinom
#' @keywords distribution
#'
#' @export

dbnbinom <- function(x, size, alpha = 1, beta = 1, log = FALSE) {
  .Call('extraDistr_cpp_dbnbinom', PACKAGE = 'extraDistr', x, size, alpha, beta, log)
}


#' @rdname BetaNegBinom
#' @export

rbnbinom <- function(n, size, alpha = 1, beta = 1) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rbnbinom', PACKAGE = 'extraDistr', n, size, alpha, beta)
}

