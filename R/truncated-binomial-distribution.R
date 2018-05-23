

#' Truncated binomial distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the truncated binomial distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param size	          number of trials (zero or more).
#' @param prob            probability of success on each trial.
#' @param a,b             lower and upper truncation points (\code{a < x <= b}).
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#' 
#' @examples 
#' 
#' x <- rtbinom(1e5, 100, 0.83, 76, 86)
#' xx <- seq(0, 100)
#' plot(prop.table(table(x)))
#' lines(xx, dtbinom(xx, 100, 0.83, 76, 86), col = "red")
#' hist(ptbinom(x, 100, 0.83, 76, 86))
#' 
#' xx <- seq(0, 100, by = 0.01)
#' plot(ecdf(x))
#' lines(xx, ptbinom(xx, 100, 0.83, 76, 86), col = "red", lwd = 2)
#' uu <- seq(0, 1, by = 0.001)
#' lines(qtbinom(uu, 100, 0.83, 76, 86), uu, col = "blue", lty = 2)
#'
#' @name TruncBinom
#' @aliases TruncBinom
#' @aliases dtbinom
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Discrete
#'
#' @export

dtbinom <- function(x, size, prob, a = -Inf, b = Inf, log = FALSE) {
  cpp_dtbinom(x, size, prob, a, b, log[1L])
}


#' @rdname TruncBinom
#' @export

ptbinom <- function(q, size, prob, a = -Inf, b = Inf, lower.tail = TRUE, log.p = FALSE) {
  cpp_ptbinom(q, size, prob, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname TruncBinom
#' @export

qtbinom <- function(p, size, prob, a = -Inf, b = Inf, lower.tail = TRUE, log.p = FALSE) {
  cpp_qtbinom(p, size, prob, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname TruncBinom
#' @export

rtbinom <- function(n, size, prob, a = -Inf, b = Inf) {
  if (length(n) > 1) n <- length(n)
  cpp_rtbinom(n, size, prob, a, b)
}

