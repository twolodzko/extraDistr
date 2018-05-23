

#' Inverse chi-squared and scaled chi-squared distributions
#'
#' Density, distribution function and random generation
#' for the inverse chi-squared distribution and scaled chi-squared distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param nu              positive valued shape parameter.
#' @param tau             positive valued scaling parameter; if provided it
#'                        returns values for scaled chi-squared distributions.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' If \eqn{X} follows \eqn{\chi^2 (\nu)} distribution, then \eqn{1/X} follows inverse
#' chi-squared distribution parametrized by \eqn{\nu}. Inverse chi-squared distribution
#' is a special case of inverse gamma distribution with parameters
#' \eqn{\alpha=\frac{\nu}{2}}{\alpha=\nu/2} and \eqn{\beta=\frac{1}{2}}{\beta=1/2};
#' or \eqn{\alpha=\frac{\nu}{2}}{\alpha=\nu/2} and
#' \eqn{\beta=\frac{\nu\tau^2}{2}}{\beta=(\nu\tau^2)/2} for scaled inverse
#' chi-squared distribution.
#'
#' @seealso \code{\link[stats]{Chisquare}}, \code{\link[stats]{GammaDist}}
#' 
#' @examples 
#' 
#' x <- rinvchisq(1e5, 20)
#' hist(x, 100, freq = FALSE)
#' curve(dinvchisq(x, 20), 0, 1, n = 501, col = "red", add = TRUE)
#' hist(pinvchisq(x, 20))
#' plot(ecdf(x))
#' curve(pinvchisq(x, 20), 0, 1, n = 501, col = "red", lwd = 2, add = TRUE)
#' 
#' # scaled
#' 
#' x <- rinvchisq(1e5, 10, 5)
#' hist(x, 100, freq = FALSE)
#' curve(dinvchisq(x, 10, 5), 0, 150, n = 501, col = "red", add = TRUE)
#' hist(pinvchisq(x, 10, 5))
#' plot(ecdf(x))
#' curve(pinvchisq(x, 10, 5), 0, 150, n = 501, col = "red", lwd = 2, add = TRUE)
#'
#' @name InvChiSq
#' @aliases InvChiSq
#' @aliases dinvchisq
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dinvchisq <- function(x, nu, tau, log = FALSE) {
  if (missing(tau))
    dinvgamma(x, nu/2, 0.5, log = log[1L])
  else
    dinvgamma(x, nu/2, (nu*tau)/2, log = log[1L])
}


#' @rdname InvChiSq
#' @export

pinvchisq <- function(q, nu, tau, lower.tail = TRUE, log.p = FALSE) {
  if (missing(tau))
    cpp_pinvgamma(q, nu/2, 0.5, lower.tail, log.p[1L])
  else
    cpp_pinvgamma(q, nu/2, (nu*tau)/2, lower.tail, log.p[1L])
}


#' @rdname InvChiSq
#' @export

qinvchisq <- function(p, nu, tau, lower.tail = TRUE, log.p = FALSE) {
  if (missing(tau))
    1/qchisq(p, nu, lower.tail = !lower.tail[1L], log.p = log.p[1L])
  else
    1/qgamma(p, nu/2, (nu*tau)/2, lower.tail = !lower.tail[1L], log.p = log.p[1L])
}


#' @rdname InvChiSq
#' @export

rinvchisq <- function(n, nu, tau) {
  if (length(n) > 1) 
    n <- length(n)
  if (missing(tau))
    1/rchisq(n, nu)
  else
    1/rgamma(n, nu/2, (nu*tau)/2)
}

