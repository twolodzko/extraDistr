

#' Bernoulli distribution
#'
#' Probability mass function, distribution function, quantile function and random generation
#' for the Bernoulli distribution.
#'
#' @param x 	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param prob            probability of success.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @name Bernoulli
#' @aliases Bernoulli
#' @aliases dbern
#' @keywords distribution
#'
#' @export

dbern <- function(x, prob = 0.5, log = FALSE) {
  .Call('extraDistr_cpp_dbern', PACKAGE = 'extraDistr', x, prob, log)
}


#' @rdname Bernoulli
#' @export

pbern <- function(x, prob = 0.5, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_pbern', PACKAGE = 'extraDistr', x, prob, lower.tail, log.p)
}


#' @rdname Bernoulli
#' @export

qbern <- function(p, prob = 0.5, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qbern', PACKAGE = 'extraDistr', p, prob, lower.tail, log.p)
}


#' @rdname Bernoulli
#' @export

rbern <- function(n, prob = 0.5) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rbern', PACKAGE = 'extraDistr', n, prob)
}

