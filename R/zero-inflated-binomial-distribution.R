

#' Zero-inflated binomial distribution
#'
#' Probability mass function and random generation
#' for the zero-inflated binomial distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param size	          number of trials (zero or more).
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
#' \pi + (1 - \pi) (1-p)^n & x = 0 \\
#' (1 - \pi) {n \choose x} p^x (1-p)^{n-x} & x > 0 \\
#' \end{array}\right.
#' }{
#' f(x) = [if x = 0:] (1-\pi)+\pi * p^r [else:] (1-\pi) * dnbinom(x, size, prob)
#' }
#' 
#' @seealso \code{\link[stats]{Binomial}}
#' 
#' @examples 
#' 
#' x <- rzib(1e5, 10, 0.6, 0.33)
#' xx <- -2:20
#' plot(prop.table(table(x)), type = "h")
#' lines(xx, dzib(xx, 10, 0.6, 0.33), col = "red")
#' 
#' xx <- seq(0, 20, by = 0.01)
#' plot(ecdf(x))
#' lines(xx, pzib(xx, 10, 0.6, 0.33), col = "red")
#'
#' @name ZIB
#' @aliases ZIB
#' @aliases dzib
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Discrete
#'
#' @export

dzib <- function(x, size, prob, pi, log = FALSE) {
  cpp_dzib(x, size, prob, pi, log[1L])
}


#' @rdname ZIB
#' @export

pzib <- function(q, size, prob, pi, lower.tail = TRUE, log.p = FALSE) {
  cpp_pzib(q, size, prob, pi, lower.tail[1L], log.p[1L])
}


#' @rdname ZIB
#' @export

qzib <- function(p, size, prob, pi, lower.tail = TRUE, log.p = FALSE) {
  cpp_qzib(p, size, prob, pi, lower.tail[1L], log.p[1L])
}


#' @rdname ZIB
#' @export

rzib <- function(n, size, prob, pi) {
  if (length(n) > 1) n <- length(n)
  cpp_rzib(n, size, prob, pi)
}

