

#' Nakagami distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Nakagami distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param m               positive valued shape parameter.
#' @param w               positive valued scale parameter.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' The Nakagami distribution (Nakagami, 1960) with shape \eqn{m} and scale
#' \eqn{\omega} has density
#' \deqn{
#' f(x) = \frac{2m^m}{\Gamma(m)\omega^m} x^{2m-1} \exp\left(-\frac{m}{\omega} x^2\right)
#' }{
#' f(x) = (2*m^m / (gamma(m)*\omega^m)) * x^(2m-1) * exp(-m/\omega * x^2)
#' }
#' for \eqn{x \ge 0}, \eqn{m > 0} and \eqn{\omega > 0}.
#'
#' If \eqn{Y} is Gamma distributed with shape \eqn{m} and rate \eqn{m/\omega},
#' then \eqn{X = \sqrt{Y}} is Nakagami distributed with shape \eqn{m} and
#' scale \eqn{\omega}.
#'
#' @references
#' Nakagami, M. (1960). The m-Distribution - A General Formula of
#' Intensity Distribution of Rapid Fading. In Statistical Methods in
#' Radio Wave Propagation, edited by W. C. Hoffman, 3-36. Pergamon Press.
#'
#' @seealso \code{\link[stats]{GammaDist}}
#'
#' @examples
#'
#' x <- rnaka(1e5, m = 2, w = 1)
#' hist(x, 100, freq = FALSE)
#' curve(dnaka(x, m = 2, w = 1), 0, 3, col = "red", add = TRUE)
#' hist(pnaka(x, m = 2, w = 1))
#' plot(ecdf(x))
#' curve(pnaka(x, m = 2, w = 1), 0, 3, col = "red", lwd = 2, add = TRUE)
#'
#' @name Nakagami
#' @aliases Nakagami
#' @aliases dnaka
#'
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dnaka <- function(x, m, w, log = FALSE) {
  cpp_dnaka(x, m, w, log[1L])
}


#' @rdname Nakagami
#' @export

pnaka <- function(q, m, w, lower.tail = TRUE, log.p = FALSE) {
  cpp_pnaka(q, m, w, lower.tail[1L], log.p[1L])
}


#' @rdname Nakagami
#' @export

qnaka <- function(p, m, w, lower.tail = TRUE, log.p = FALSE) {
  cpp_qnaka(p, m, w, lower.tail[1L], log.p[1L])
}


#' @rdname Nakagami
#' @export

rnaka <- function(n, m, w) {
  if (length(n) > 1) n <- length(n)
  cpp_rnaka(n, m, w)
}
