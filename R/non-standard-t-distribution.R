

#' Non-standard t-distribution
#' 
#' Probability mass function, distribution function and random generation
#' for non-standard t-distribution. Non-standard t-distribution besides
#' degrees of freedom \eqn{\nu}, is parametrized using additional parameters
#' \eqn{\mu} for location and \eqn{\sigma} for scale (\eqn{\mu=0} and
#' \eqn{\sigma = 1} for standard t-distribution).
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param mu              vector of locations
#' @param sigma           vector of positive valued scale parameters.
#' @param df              degrees of freedom (> 0, maybe non-integer). \code{df = Inf} is allowed.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#' 
#' @seealso \code{\link[stats]{TDist}}
#' 
#' @examples 
#' 
#' x <- rnst(1e5, 1000, 5, 13)
#' hist(x, 100, freq = FALSE)
#' curve(dnst(x, 1000, 5, 13), -60, 60, col = "red", add = TRUE)
#' hist(pnst(x, 1000, 5, 13))
#' plot(ecdf(x))
#' curve(pnst(x, 1000, 5, 13), -60, 60, col = "red", lwd = 2, add = TRUE)
#' 
#' @name NonStandardT
#' @aliases NonStandardT
#' @aliases dnst
#' @keywords distribution
#' 
#' @export

dnst <- function(x, df, mu = 0, sigma = 1, log = FALSE) {
  cpp_dnst(x, df, mu, sigma, log)
}


#' @rdname NonStandardT
#' @export

pnst <- function(q, df, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_pnst(q, df, mu, sigma, lower.tail, log.p)
}


#' @rdname NonStandardT
#' @export

qnst <- function(p, df, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_qnst(p, df, mu, sigma, lower.tail, log.p)
}


#' @rdname NonStandardT
#' @export

rnst <- function(n, df, mu = 0, sigma = 1) {
  if (length(n) > 1) n <- length(n)
  cpp_rnst(n, df, mu, sigma)
}

