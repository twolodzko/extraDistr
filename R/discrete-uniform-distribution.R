

#' Discrete uniform distribution
#'
#' Probability mass function, distribution function, quantile function and random generation
#' for the Bernoulli distribution.
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
#' @examples 
#' 
#' x <- sample.int(10, 1e5, replace = TRUE) 
#' xx <- -1:11
#' plot(prop.table(table(x)), type = "h")
#' lines(xx, ddunif(xx, 1, 10), col = "red")
#' hist(pdunif(x, 1, 10))
#' plot(ecdf(x))
#' lines(xx, pdunif(xx, 1, 10), col = "red")
#'
#' @name DiscreteUniform
#' @aliases DiscreteUniform
#' @aliases ddunif
#' @keywords distribution
#'
#' @export

ddunif <- function(x, min, max, log = FALSE) {
  .Call('extraDistr_cpp_ddunif', PACKAGE = 'extraDistr', x, min, max, log)
}


#' @rdname DiscreteUniform
#' @export

pdunif <- function(q, min, max, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_pdunif', PACKAGE = 'extraDistr', q, min, max, lower.tail, log.p)
}


#' @rdname DiscreteUniform
#' @export

qdunif <- function(p, min, max, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qdunif', PACKAGE = 'extraDistr', p, min, max, lower.tail, log.p)
}


#' @rdname DiscreteUniform
#' @export

rdunif <- function(n, min, max) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rdunif', PACKAGE = 'extraDistr', n, min, max)
}

