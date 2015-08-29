

#' Dirichlet distribution
#'
#' Density function, cumulative distribution function and random generation
#' for the Dirichlet distribution.
#'
#' @param x               matrix of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param alpha           vector; concentration parameter.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#' @param nsim            number of samples in Monte Carlo simulation for calculating
#'                        cumulative distribution function; the higher is more precise
#'                        but slower.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = \frac{\Gamma(\sum_k \alpha_k)}{\prod_k \Gamma(\alpha_k)} \prod_k x_k^{k-1}
#' }{
#' f(x) = \Gamma(sum(\alpha[k])) / prod(\Gamma(\alpha[k])) * prod(x[k]^{k-1})
#' }
#'
#' Cumulative distribution function is approximated using Monte Carlo simulation.
#'
#' @references
#' Devroye, L. (1986). Non-Uniform Random Variate Generation. Springer-Verlag.
#'
#' @references
#' \url{http://stats.stackexchange.com/a/57313/35989}
#'
#' @name Dirichlet
#' @aliases Dirichlet
#' @aliases ddirichlet
#' @export

ddirichlet <- function(x, alpha, log = FALSE) {
  if (is.data.frame(x)) x <- as.matrix(x)
  .Call('extraDistr_cpp_ddirichlet', PACKAGE = 'extraDistr', x, alpha, log)
}


#' @rdname Dirichlet
#' @export

pdirichlet <- function(x, alpha, lower.tail = TRUE, log.p = FALSE, nsim = 10000L) {
  if (is.data.frame(x)) x <- as.matrix(x)
  .Call('extraDistr_cpp_pdirichlet', PACKAGE = 'extraDistr', x, alpha, lower.tail, log.p, nsim)
}


#' @rdname Dirichlet
#' @export

rdirichlet <- function (n, alpha) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rdirichlet', PACKAGE = 'extraDistr', n, alpha)
}

