

#' Slash distribution
#' 
#' Probability mass function, distribution function and random generation
#' for slash distribution.
#'
#' @param x 	            vector of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param mu              vector of locations
#' @param sigma           vector of scale parameters.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#' 
#' @details
#' 
#' If \eqn{Z \sim \mathrm{Normal}(0, 1)}{Z ~ Normal(0, 1)} and \eqn{U \sim \mathrm{Uniform}(0, 1)}{U ~ Uniform(0, 1)},
#' then \eqn{Z/U} follows slash distribution.
#' 
#' Probability density function
#' 
#' \deqn{
#' f(x) = \frac{\phi(0) - \phi(x)}{x^2}
#' }{
#' f(x) = [\phi(0)-\phi(x)]/x^2
#' }
#' 
#' Cumulative distribution function
#' 
#' \deqn{
#' F(x) = \left\{\begin{array}{ll}
#' \Phi(x) - \frac{\phi(0)-\phi(x)}{x} & x \neq 0 \\
#' \frac{1}{2} & x = 0
#' \end{array}\right.
#' }{
#' F(x) = [if x != 0:] \Phi(x) - [\phi(0)-\phi(x)]/x [else:] 1/2
#' }
#' 
#' @name Slash
#' @aliases Slash
#' @aliases dslash
#' @keywords distribution
#' 
#' @export

dslash <- function(x, mu = 0, sigma = 1, log = FALSE) {
  .Call('extraDistr_cpp_dslash', PACKAGE = 'extraDistr', x, mu, sigma, log)
}


#' @rdname Slash
#' @export

pslash <- function(x, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_pslash', PACKAGE = 'extraDistr', x, mu, sigma, lower.tail, log.p)
}


#' @rdname Slash
#' @export

rslash <- function(n, mu = 0, sigma = 1) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rslash', PACKAGE = 'extraDistr', n, mu, sigma)
}

