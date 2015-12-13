

#' Inverse-chi-squared distribution
#'
#' Density, distribution function and random generation
#' for the Inverse-chi-squared distribution.
#'
#' @param x 	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param nu,tau          shape and scale parameters.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability mass function
#' \deqn{
#' f(x) = \frac{2^{-\nu/2}}{\Gamma(\nu/2)} x^{-\nu/2-1} e^{-1/2x}
#' }{
#' f(x) = (2^(-\nu/2))/\Gamma(\nu/2) * x^(-\nu/2-1) * exp(-1/(2*x))
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = \frac{\gamma(\nu/2, 1/2x)}{\Gamma(\nu/2)}
#' }{
#' F(x) = \gamma(\nu/2, 1/(2*x)) / \Gamma(\nu/2)
#' }
#'
#' @seealso \code{\link{Chisquare}}
#'
#' @name ScInvChiSq
#' @aliases ScInvChiSq
#' @aliases dsinvchisq
#' @keywords distribution
#'
#' @export

dsinvchisq <- function(x, nu, tau = 1, log = FALSE) {
  .Call('extraDistr_cpp_dsinvchisq', PACKAGE = 'extraDistr', x, nu, tau, log)
}


#' @rdname ScInvChiSq
#' @export

psinvchisq <- function(x, nu, tau = 1, lower.tail = TRUE, log.p = FALSE) {
  pgamma((tau^2 * nu)/(2*x), nu/2, 1, lower.tail = !lower.tail, log.p = log.p)
}


#' @rdname ScInvChiSq
#' @export

rsinvchisq <- function(n, nu, tau = 1) {
  if (length(n) > 1) n <- length(n)
  1/rgamma(n, nu/2, (nu*tau^2)/2)
}

