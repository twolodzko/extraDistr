

#' Discrete uniform distribution
#'
#' Probability mass function, distribution function, quantile function and random generation
#' for the discrete uniform distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param min,max         lower and upper limits of the distribution. Must be finite.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#' 
#' @details 
#' 
#' If \code{min == max}, then discrete uniform distribution is a degenerate distribution.
#'                                           
#' @examples 
#' 
#' x <- rdunif(1e5, 1, 10) 
#' xx <- -1:11
#' plot(prop.table(table(x)), type = "h")
#' lines(xx, ddunif(xx, 1, 10), col = "red")
#' hist(pdunif(x, 1, 10))
#' xx <- seq(-1, 11, by = 0.01)
#' plot(ecdf(x))
#' lines(xx, pdunif(xx, 1, 10), col = "red")
#'
#' @name DiscreteUniform
#' @aliases DiscreteUniform
#' @aliases ddunif
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Discrete
#'
#' @export

ddunif <- function(x, min, max, log = FALSE) {
  cpp_ddunif(x, min, max, log[1L])
}


#' @rdname DiscreteUniform
#' @export

pdunif <- function(q, min, max, lower.tail = TRUE, log.p = FALSE) {
  cpp_pdunif(q, min, max, lower.tail[1L], log.p[1L])
}


#' @rdname DiscreteUniform
#' @export

qdunif <- function(p, min, max, lower.tail = TRUE, log.p = FALSE) {
  cpp_qdunif(p, min, max, lower.tail[1L], log.p[1L])
}


#' @rdname DiscreteUniform
#' @export

rdunif <- function(n, min, max) {
  if (length(n) > 1) n <- length(n)
  cpp_rdunif(n, min, max)
}

