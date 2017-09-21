

#' Location-scale version of the t-distribution
#' 
#' Probability mass function, distribution function and random generation
#' for location-scale version of the t-distribution. Location-scale version
#' of the t-distribution besides degrees of freedom \eqn{\nu}, is parametrized
#' using additional parameters \eqn{\mu} for location and \eqn{\sigma} for
#' scale (\eqn{\mu=0} and \eqn{\sigma = 1} for standard t-distribution).
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
#' x <- rlst(1e5, 1000, 5, 13)
#' hist(x, 100, freq = FALSE)
#' curve(dlst(x, 1000, 5, 13), -60, 60, col = "red", add = TRUE)
#' hist(plst(x, 1000, 5, 13))
#' plot(ecdf(x))
#' curve(plst(x, 1000, 5, 13), -60, 60, col = "red", lwd = 2, add = TRUE)
#' 
#' @name LocationScaleT
#' @aliases LocationScaleT
#' @aliases dlst
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#' 
#' @export

dlst <- function(x, df, mu = 0, sigma = 1, log = FALSE) {
  cpp_dlst(x, df, mu, sigma, log[1L])
}


#' @rdname LocationScaleT
#' @export

plst <- function(q, df, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_plst(q, df, mu, sigma, lower.tail[1L], log.p[1L])
}


#' @rdname LocationScaleT
#' @export

qlst <- function(p, df, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_qlst(p, df, mu, sigma, lower.tail[1L], log.p[1L])
}


#' @rdname LocationScaleT
#' @export

rlst <- function(n, df, mu = 0, sigma = 1) {
  if (length(n) > 1) n <- length(n)
  cpp_rlst(n, df, mu, sigma)
}





# Depreciated since ver 1.8.7
# To be removed in ver 1.8.8

#' Deprecated functions
#' 
#' \emph{Warning}: the \code{dnst}, \code{pnst}, \code{qnst}, and \code{rnst}
#' functions are now depreciated, please use \code{dlst}, \code{plst}, \code{qlst},
#' and \code{rlst} functions instead.
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
#' @name extraDistr-deprecated
#' @aliases dnst
#' 
#' @export

dnst <- function(x, df, mu = 0, sigma = 1, log = FALSE) {
  .Deprecated("dlst")
  dlst(x, df, mu, sigma, log[1L])
}

#' @rdname extraDistr-deprecated
#' @export

pnst <- function(q, df, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .Deprecated("plst")
  plst(q, df, mu, sigma, lower.tail[1L], log.p[1L])
}

#' @rdname extraDistr-deprecated
#' @export

qnst <- function(p, df, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .Deprecated("qlst")
  qlst(p, df, mu, sigma, lower.tail[1L], log.p[1L])
}

#' @rdname extraDistr-deprecated
#' @export

rnst <- function(n, df, mu = 0, sigma = 1) {
  .Deprecated("rlst")
  rlst(n, df, mu, sigma)
}

