

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
#' f(x) = \frac{x^{-\alpha-1} \exp(-\frac{1}{\beta x})}{\Gamma(\alpha) \beta^\alpha}
#' }{
#' f(x) = (x^(-\alpha-1) * exp(-1/(\beta*x))) / (\Gamma(\alpha)*\beta^\alpha)
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = \frac{\gamma(\alpha, \frac{1}{\beta x})}{\Gamma(\alpha)}
#' }{
#' F(x) = \gamma(\alpha, 1/(\beta*x)) / \Gamma(\alpha)
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
#' xx <- seq(0, 1, by = 0.001)
#' hist(x, 100, freq = FALSE)
#' lines(xx, dinvgamma(xx, 20, 3), col = "red")
#' hist(pinvgamma(x, 20, 3))
#' plot(ecdf(x))
#' lines(xx, pinvgamma(xx, 20, 3), col = "red", lwd = 2) 
#'
#' @name InvGamma
#' @aliases InvGamma
#' @aliases dinvgamma
#' @keywords distribution
#'
#' @export

dinvgamma <- function(x, alpha, beta = 1, log = FALSE) {
  cpp_dinvgamma(x, alpha, 1/beta, log)
}


#' @rdname InvGamma
#' @export

pinvgamma <- function(q, alpha, beta = 1, lower.tail = TRUE, log.p = FALSE) {
  pgamma(1/q, alpha, beta, lower.tail = !lower.tail, log.p = log.p)
}


#' @rdname InvGamma
#' @export

qinvgamma <- function(p, alpha, beta = 1, lower.tail = TRUE, log.p = FALSE) {
  1/qgamma(p, alpha, beta, lower.tail = !lower.tail, log.p = log.p)
}


#' @rdname InvGamma
#' @export

rinvgamma <- function(n, alpha, beta = 1) {
  1/rgamma(n, alpha, beta)
}

