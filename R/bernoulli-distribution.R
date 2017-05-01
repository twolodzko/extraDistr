

#' Bernoulli distribution
#'
#' Probability mass function, distribution function, quantile function and random generation
#' for the Bernoulli distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param prob            probability of success; (\code{0 < prob < 1}).
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'                        
#' @seealso \code{\link[stats]{Binomial}}
#' 
#' @examples 
#' 
#' prop.table(table(rbern(1e5, 0.5)))
#'
#' @name Bernoulli
#' @aliases Bernoulli
#' @aliases dbern
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Discrete
#'
#' @export

dbern <- function(x, prob = 0.5, log = FALSE) {
  cpp_dbern(x, prob, log[1L])
}


#' @rdname Bernoulli
#' @export

pbern <- function(q, prob = 0.5, lower.tail = TRUE, log.p = FALSE) {
  cpp_pbern(q, prob, lower.tail[1L], log.p[1L])
}


#' @rdname Bernoulli
#' @export

qbern <- function(p, prob = 0.5, lower.tail = TRUE, log.p = FALSE) {
  cpp_qbern(p, prob, lower.tail[1L], log.p[1L])
}


#' @rdname Bernoulli
#' @export

rbern <- function(n, prob = 0.5) {
  if (length(n) > 1) n <- length(n)
  cpp_rbern(n, prob)
}

