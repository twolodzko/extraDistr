

#' Zero-inflated Poisson distribution
#'
#' Probability mass function and random generation
#' for the zero-inflated Poisson distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param lambda	        vector of (non-negative) means.
#' @param pi              probability of extra zeros.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = \left\{\begin{array}{ll}
#' \pi + (1 - \pi) e^{-\lambda} & x = 0 \\
#' (1 - \pi) \frac{\lambda^{x} e^{-\lambda}} {x!} & x > 0 \\
#' \end{array}\right.
#' }{
#' f(x) = [if x = 0:] \pi + (1-\pi) * exp(-\lambda) [else:] (1-\pi) * dpois(x, lambda)
#' }
#' 
#' @seealso \code{\link[stats]{Poisson}}
#' 
#' @examples 
#' 
#' x <- rzip(1e5, 6, 0.33)
#' xx <- -2:20
#' plot(prop.table(table(x)), type = "h")
#' lines(xx, dzip(xx, 6, 0.33), col = "red")
#' 
#' xx <- seq(0, 20, by = 0.01)
#' plot(ecdf(x))
#' lines(xx, pzip(xx, 6, 0.33), col = "red")
#'
#' @name ZIP
#' @aliases ZIP
#' @aliases dzip
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Discrete
#'
#' @export

dzip <- function(x, lambda, pi, log = FALSE) {
  cpp_dzip(x, lambda, pi, log[1L])
}


#' @rdname ZIP
#' @export

pzip <- function(q, lambda, pi, lower.tail = TRUE, log.p = FALSE) {
  cpp_pzip(q, lambda, pi, lower.tail[1L], log.p[1L])
}


#' @rdname ZIP
#' @export

qzip <- function(p, lambda, pi, lower.tail = TRUE, log.p = FALSE) {
  cpp_qzip(p, lambda, pi, lower.tail[1L], log.p[1L])
}


#' @rdname ZIP
#' @export

rzip <- function(n, lambda, pi) {
  if (length(n) > 1) n <- length(n)
  cpp_rzip(n, lambda, pi)
}

