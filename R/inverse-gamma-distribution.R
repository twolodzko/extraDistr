

#' Inverse-gamma distribution
#'
#' Density, distribution function and random generation
#' for the inverse-gamma distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param alpha,beta      positive valued shape and scale parameters.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability mass function
#' \deqn{
#' f(x) = \frac{\beta^\alpha x^{-\alpha-1} \exp(-\frac{\beta}{x})}{\Gamma(\alpha)}
#' }{
#' f(x) = (\beta^\alpha * x^(-\alpha-1) * exp(-\beta/x)) / \Gamma(\alpha)
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = \frac{\gamma(\alpha, \frac{\beta}{x})}{\Gamma(\alpha)}
#' }{
#' F(x) = \gamma(\alpha, \beta/x) / \Gamma(\alpha)
#' }
#'
#' @references
#' Witkovsky, V. (2001). Computing the distribution of a linear
#' combination of inverted gamma variables. Kybernetika 37(1), 79-90.
#'
#' @references
#' Leemis, L.M. and McQueston, L.T. (2008). Univariate Distribution
#' Relationships. American Statistician 62(1): 45-53.
#'
#' @seealso \code{\link[stats]{GammaDist}}
#' 
#' @examples 
#' 
#' x <- rinvgamma(1e5, 20, 3)
#' hist(x, 100, freq = FALSE)
#' curve(dinvgamma(x, 20, 3), 0, 1, col = "red", add = TRUE, n = 5000)
#' hist(pinvgamma(x, 20, 3))
#' plot(ecdf(x))
#' curve(pinvgamma(x, 20, 3), 0, 1, col = "red", lwd = 2, add = TRUE, n = 5000) 
#'
#' @name InvGamma
#' @aliases InvGamma
#' @aliases dinvgamma
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dinvgamma <- function(x, alpha, beta = 1, log = FALSE) {
  cpp_dinvgamma(x, alpha, 1/beta, log[1L])
}


#' @rdname InvGamma
#' @export

pinvgamma <- function(q, alpha, beta = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_pinvgamma(q, alpha, beta, lower.tail, log.p[1L])
}


#' @rdname InvGamma
#' @export

qinvgamma <- function(p, alpha, beta = 1, lower.tail = TRUE, log.p = FALSE) {
  1/qgamma(p, alpha, beta, lower.tail = !lower.tail[1L], log.p = log.p[1L])
}


#' @rdname InvGamma
#' @export

rinvgamma <- function(n, alpha, beta = 1) {
  1/rgamma(n, alpha, beta)
}

