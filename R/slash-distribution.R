

#' Slash distribution
#' 
#' Probability mass function, distribution function and random generation
#' for slash distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param mu              vector of locations
#' @param sigma           vector of positive valued scale parameters.
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
#' f(x) = \left\{\begin{array}{ll}
#' \frac{\phi(0) - \phi(x)}{x^2} & x \ne 0 \\
#' \frac{1}{2\sqrt{2\pi}} & x = 0
#' \end{array}\right.
#' }{
#' f(x) = [if x != 0:] (\phi(0)-\phi(x))/x^2 [else:] 1/(2*sqrt(2*\pi))
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
#' @examples 
#' 
#' x <- rslash(1e5, 5, 3)
#' hist(x, 1e5, freq = FALSE, xlim = c(-100, 100))
#' curve(dslash(x, 5, 3), -100, 100, col = "red", n = 500, add = TRUE)
#' hist(pslash(x, 5, 3))
#' plot(ecdf(x), xlim = c(-100, 100))
#' curve(pslash(x, 5, 3), -100, 100, col = "red", lwd = 2, n = 500, add = TRUE)
#' 
#' @name Slash
#' @aliases Slash
#' @aliases dslash
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#' 
#' @export

dslash <- function(x, mu = 0, sigma = 1, log = FALSE) {
  cpp_dslash(x, mu, sigma, log[1L])
}


#' @rdname Slash
#' @export

pslash <- function(q, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_pslash(q, mu, sigma, lower.tail[1L], log.p[1L])
}


#' @rdname Slash
#' @export

rslash <- function(n, mu = 0, sigma = 1) {
  if (length(n) > 1) n <- length(n)
  cpp_rslash(n, mu, sigma)
}

