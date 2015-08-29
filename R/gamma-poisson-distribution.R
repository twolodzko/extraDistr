

#' Gamma-Poisson distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Gamma-Poisson distribution.
#'
#' @param x 	            vector of quantiles.
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
#' f(x) = \frac{\Gamma(x+\beta) \alpha^x}{\Gamma(\beta) (1+\alpha)^{\beta+x} x!}
#' }{
#' f(x) = (\Gamma(x+\beta)*\alpha^x) / (\Gamma(\beta)*(1+\alpha)^(\beta+x) * x!)
#' }
#'
#' @name GammaPoiss
#' @aliases GammaPoiss
#' @aliases dgpois
#' @keywords distribution
#'
#' @export

dgpois <- function(x, alpha, beta, log = FALSE) {
  .Call('extraDistr_cpp_dgpois', PACKAGE = 'extraDistr', x, alpha, beta, log)
}


#' @rdname GammaPoiss
#' @export

pgpois <- function(x, alpha, beta, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_pgpois', PACKAGE = 'extraDistr', x, alpha, beta, lower.tail, log.p)
}


#' @rdname GammaPoiss
#' @export

rgpois <- function(n, alpha, beta) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rgpois', PACKAGE = 'extraDistr', n, alpha, beta)
}

