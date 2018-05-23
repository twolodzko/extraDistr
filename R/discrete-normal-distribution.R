

#' Discrete normal distribution
#' 
#' Probability mass function, distribution function and random generation
#' for discrete normal distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param mean            vector of means.
#' @param sd              vector of standard deviations.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#' @details
#' 
#' Probability mass function
#' 
#' \deqn{
#' f(x) = \Phi\left(\frac{x-\mu+1}{\sigma}\right) - \Phi\left(\frac{x-\mu}{\sigma}\right)
#' }{
#' f(x) = \Phi((x-\mu+1)/\sigma) - \Phi((x-\mu)/\sigma)
#' }
#' 
#' Cumulative distribution function
#' 
#' \deqn{
#' F(x) = \Phi\left(\frac{\lfloor x \rfloor + 1 - \mu}{\sigma}\right)
#' }{
#' F(x) = \Phi((floor(x)+1-\mu)/\sigma)
#' }
#' 
#' @references 
#' Roy, D. (2003). The discrete normal distribution.
#' Communications in Statistics-Theory and Methods, 32, 1871-1883.
#' 
#' @seealso \code{\link[stats]{Normal}}
#' 
#' @examples 
#' 
#' x <- rdnorm(1e5, 0, 3)
#' xx <- -15:15
#' plot(prop.table(table(x)))
#' lines(xx, ddnorm(xx, 0, 3), col = "red")
#' hist(pdnorm(x, 0, 3))
#' plot(ecdf(x))
#' xx <- seq(-15, 15, 0.1)
#' lines(xx, pdnorm(xx, 0, 3), col = "red", lwd = 2, type = "s")
#' 
#' @name DiscreteNormal
#' @aliases DiscreteNormal
#' @aliases ddnorm
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Discrete
#' 
#' @export

ddnorm <- function(x, mean = 0, sd = 1, log = FALSE) {
  cpp_ddnorm(x, mean, sd, log[1L])
}


#' @rdname DiscreteNormal
#' @export

pdnorm <- function(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) {
  pnorm(floor(q)+1, mean, sd, lower.tail[1L], log.p[1L])
}


#' @rdname DiscreteNormal
#' @export

rdnorm <- function(n, mean = 0, sd = 1) {
  floor(rnorm(n, mean, sd))
}

