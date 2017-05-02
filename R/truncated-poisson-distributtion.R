

#' Truncated Poisson distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the truncated Poisson distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param lambda	        vector of (non-negative) means.
#' @param a,b             lower and upper truncation points (\code{a < x <= b}).
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @references
#' Plackett, R.L. (1953). The truncated Poisson distribution.
#' Biometrics, 9(4), 485-488.
#' 
#' @references 
#' Singh, J. (1978). A characterization of positive Poisson distribution and
#' its statistical application.
#' SIAM Journal on Applied Mathematics, 34(3), 545-548.
#' 
#' @references 
#' Dalgaard, P. (May 1, 2005). [R] simulate zero-truncated Poisson distribution.
#' R-help mailing list.
#' \url{https://stat.ethz.ch/pipermail/r-help/2005-May/070680.html}
#' 
#' @examples 
#' 
#' x <- rtpois(1e5, 14, 16)
#' xx <- seq(-1, 50)
#' plot(prop.table(table(x)))
#' lines(xx, dtpois(xx, 14, 16), col = "red")
#' hist(ptpois(x, 14, 16))
#' 
#' xx <- seq(0, 50, by = 0.01)
#' plot(ecdf(x))
#' lines(xx, ptpois(xx, 14, 16), col = "red", lwd = 2)
#' 
#' uu <- seq(0, 1, by = 0.001)
#' lines(qtpois(uu, 14, 16), uu, col = "blue", lty = 2)
#' 
#' # Zero-truncated Poisson
#' 
#' x <- rtpois(1e5, 5, 0)
#' xx <- seq(-1, 50)
#' plot(prop.table(table(x)))
#' lines(xx, dtpois(xx, 5, 0), col = "red")
#' hist(ptpois(x, 5, 0))
#' 
#' xx <- seq(0, 50, by = 0.01)
#' plot(ecdf(x))
#' lines(xx, ptpois(xx, 5, 0), col = "red", lwd = 2)
#' lines(qtpois(uu, 5, 0), uu, col = "blue", lty = 2)
#'
#' @name TruncPoisson
#' @aliases TruncPoisson
#' @aliases dtpois
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Discrete
#'
#' @export

dtpois <- function(x, lambda, a = -Inf, b = Inf, log = FALSE) {
  cpp_dtpois(x, lambda, a, b, log)
}


#' @rdname TruncPoisson
#' @export

ptpois <- function(q, lambda, a = -Inf, b = Inf, lower.tail = TRUE, log.p = FALSE) {
  cpp_ptpois(q, lambda, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname TruncPoisson
#' @export

qtpois <- function(p, lambda, a = -Inf, b = Inf, lower.tail = TRUE, log.p = FALSE) {
  cpp_qtpois(p, lambda, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname TruncPoisson
#' @export

rtpois <- function(n, lambda, a = -Inf, b = Inf) {
  if (length(n) > 1) n <- length(n)
  cpp_rtpois(n, lambda, a, b)
}

