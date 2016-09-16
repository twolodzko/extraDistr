

#' Discrete normal distribution
#' 
#' Probability mass function, distribution function and random generation
#' for discrete normal distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
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
#' @references 
#' Roy, D. (2003). The discrete normal distribution.
#' Communications in Statistics-Theory and Methods, 32, 1871-1883.
#' 
#' @seealso \code{\link[stats]{Normal}}
#' 
#' @examples 
#' 
#' x <- rdnorm(1e5, 7, 35)
#' xx <- -150:150
#' hist(x, 100, freq = FALSE)
#' lines(xx-0.5, ddnorm(xx, 7, 35), col = "red")
#' hist(pdnorm(x, 7, 35))
#' plot(ecdf(x))
#' lines(xx, pdnorm(xx, 7, 35), col = "red", lwd = 2)
#' 
#' @name DiscreteNormal
#' @aliases DiscreteNormal
#' @aliases ddnorm
#' @keywords distribution
#' 
#' @export

ddnorm <- function(x, mean = 0, sd = 1, log = FALSE) {
  cpp_ddnorm(x, mean, sd, log)
}


#' @rdname DiscreteNormal
#' @export

pdnorm <- function(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) {
  pnorm(q, mean, sd, lower.tail, log.p)
}


#' @rdname DiscreteNormal
#' @export

qdnorm <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) {
  ceiling(qnorm(p, mean, sd, lower.tail, log.p))
}


#' @rdname DiscreteNormal
#' @export

rdnorm <- function(n, mean = 0, sd = 1) {
  ceiling(rnorm(n, mean, sd))
}

