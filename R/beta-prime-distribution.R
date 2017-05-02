

#' Beta prime distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the beta prime distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param shape1,shape2	  non-negative parameters.
#' @param scale           positive valued scale parameter.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details 
#' 
#' If \eqn{X \sim \mathrm{Beta}(\alpha, \beta)}{X ~ Beta(\alpha, \beta)}, then
#' \eqn{\frac{X}{1-X} \sim \mathrm{BetaPrime}(\alpha, \beta)}{X/(1-X) ~ BetaPrime(\alpha, \beta)}.
#' 
#' Probability density function
#' 
#' \deqn{
#' f(x) = \frac{(x/\sigma)^{\alpha-1} (1+x/\sigma)^{-\alpha -\beta}}{\mathrm{B}(\alpha,\beta)\sigma}
#' }{
#' f(x) = ((x/\sigma)^(\alpha-1) * (1 + x/\sigma)^(-\alpha-\beta)) / (B(\alpha,\beta) * \sigma)
#' }
#' 
#' Cumulative distribution function
#' 
#' \deqn{
#' F(x) = I_{\frac{x/\sigma}{1+x/\sigma}}(\alpha, \beta)
#' }{
#' F(x) = pbeta((x/\sigma)/(1+(x/\sigma)), \alpha, \beta)
#' }
#' 
#' @seealso \code{\link[stats]{Beta}}
#'
#' @examples 
#' 
#' x <- rbetapr(1e5, 5, 3, 2)
#' hist(x, 350, freq = FALSE, xlim = c(0, 100))
#' curve(dbetapr(x, 5, 3, 2), 0, 100, col = "red", add = TRUE, n = 500)
#' hist(pbetapr(x, 5, 3, 2))
#' plot(ecdf(x), xlim = c(0, 100))
#' curve(pbetapr(x, 5, 3, 2), 0, 100, col = "red", add = TRUE, n = 500)
#' 
#' @name BetaPrime
#' @aliases BetaPrime
#' @aliases dbetapr
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dbetapr <- function(x, shape1, shape2, scale = 1, log = FALSE) {
  cpp_dbetapr(x, shape1, shape2, scale, log[1L])
}


#' @rdname BetaPrime
#' @export

pbetapr <- function(q, shape1, shape2, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_pbetapr(q, shape1, shape2, scale, lower.tail[1L], log.p[1L])
}


#' @rdname BetaPrime
#' @export

qbetapr <- function(p, shape1, shape2, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_qbetapr(p, shape1, shape2, scale, lower.tail[1L], log.p[1L])
}


#' @rdname BetaPrime
#' @export

rbetapr <- function(n, shape1, shape2, scale = 1) {
  if (length(n) > 1) n <- length(n)
  cpp_rbetapr(n, shape1, shape2, scale)
}

