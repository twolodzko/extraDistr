

#' Zero-inflated negative binomial distribution
#'
#' Probability mass function and random generation
#' for the zero-inflated negative binomial distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param size	          target for number of successful trials, or dispersion
#'                        parameter (the shape parameter of the gamma mixing
#'                        distribution). Must be strictly positive, need not be
#'                        integer.
#' @param prob            probability of success in each trial. \code{0 < prob <= 1}.
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
#' \pi + (1 - \pi) p^r & x = 0 \\
#' (1 - \pi) {x+r-1 \choose x} p^r (1-p)^x & x > 0 \\
#' \end{array}\right.
#' }{
#' f(x) = [if x = 0:] (1-\pi)+\pi * p^r [else:] (1-\pi) * dnbinom(x, size, prob)
#' }
#' 
#' @seealso \code{\link[stats]{NegBinomial}}
#' 
#' @examples 
#' 
#' x <- rzinb(1e5, 100, 0.6, 0.33)
#' xx <- -2:200
#' plot(prop.table(table(x)), type = "h")
#' lines(xx, dzinb(xx, 100, 0.6, 0.33), col = "red")
#' 
#' xx <- seq(0, 200, by = 0.01)
#' plot(ecdf(x))
#' lines(xx, pzinb(xx, 100, 0.6, 0.33), col = "red")
#'
#' @name ZINB
#' @aliases ZINB
#' @aliases dzinb
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Discrete
#'
#' @export

dzinb <- function(x, size, prob, pi, log = FALSE) {
  cpp_dzinb(x, size, prob, pi, log[1L])
}


#' @rdname ZINB
#' @export

pzinb <- function(q, size, prob, pi, lower.tail = TRUE, log.p = FALSE) {
  cpp_pzinb(q, size, prob, pi, lower.tail[1L], log.p[1L])
}


#' @rdname ZINB
#' @export

qzinb <- function(p, size, prob, pi, lower.tail = TRUE, log.p = FALSE) {
  cpp_qzinb(p, size, prob, pi, lower.tail[1L], log.p[1L])
}


#' @rdname ZINB
#' @export

rzinb <- function(n, size, prob, pi) {
  if (length(n) > 1) n <- length(n)
  cpp_rzinb(n, size, prob, pi)
}

