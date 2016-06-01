

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
#' @param s               truncation point (non-negtive); \code{s=0} (default)
#'                        for zero-truncated Poisson, otherwise values greater
#'                        than \code{s} are truncated.
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
#' xx <- seq(-1, 20)
#' plot(prop.table(table(x)))
#' lines(xx, dtpois(xx, 14, 16), col = "red")
#' hist(ptpois(x, 14, 16))
#' plot(ecdf(x))
#' lines(xx, ptpois(xx, 14, 16), col = "red", lwd = 2)
#' uu <- seq(0, 1, by = 0.001)
#' lines(qtpois(uu, 14, 16), uu, col = "blue")
#' 
#' # Zero-truncated Poisson
#' 
#' x <- rtpois(1e5, 5, 0)
#' xx <- seq(-1, 50)
#' plot(prop.table(table(x)))
#' lines(xx, dtpois(xx, 5, 0), col = "red")
#' hist(ptpois(x, 5, 0))
#' plot(ecdf(x))
#' lines(xx, ptpois(xx, 5, 0), col = "red", lwd = 2)
#' lines(qtpois(uu, 5, 0), uu, col = "blue")
#'
#' @name TruncPoisson
#' @aliases TruncPoisson
#' @aliases dtpois
#' @keywords distribution
#'
#' @export

dtpois <- function(x, lambda, s = 0, log = FALSE) {
  .Call('extraDistr_cpp_dtpois', PACKAGE = 'extraDistr', x, lambda, s, log)
}


#' @rdname TruncPoisson
#' @export

ptpois <- function(q, lambda, s = 0, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_ptpois', PACKAGE = 'extraDistr', q, lambda, s, lower.tail, log.p)
}


#' @rdname TruncPoisson
#' @export

qtpois <- function(p, lambda, s = 0, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qtpois', PACKAGE = 'extraDistr', p, lambda, s, lower.tail, log.p)
}


#' @rdname TruncPoisson
#' @export

rtpois <- function(n, lambda, s = 0) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rtpois', PACKAGE = 'extraDistr', n, lambda, s)
}

