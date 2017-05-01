

#' Bhattacharjee distribution
#'
#' Density, distribution function, and random generation for the Bhattacharjee
#' distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param mu,sigma,a	    location, scale and shape parameters.
#'                        Scale and shape must be positive.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#' 
#' If \eqn{Z \sim \mathrm{Normal}(0, 1)}{Z ~ Normal(0, 1)} and
#' \eqn{U \sim \mathrm{Uniform}(0, 1)}{U ~ Uniform(0, 1)}, then
#' \eqn{Z+U} follows Bhattacharjee distribution.
#' 
#' Probability density function
#' 
#' \deqn{
#' f(z) = \frac{1}{2a} \left[\Phi\left(\frac{x-\mu+a}{\sigma}\right) - \Phi\left(\frac{x-\mu-a}{\sigma}\right)\right]
#' }{
#' f(z) = 1/(2*a) * (\Phi((x-\mu+a)/\sigma) - \Phi((x-\mu+a)/\sigma))
#' }
#' 
#' Cumulative distribution function
#' 
#' \deqn{
#' F(z) = \frac{\sigma}{2a} \left[(x-\mu)\Phi\left(\frac{x-\mu+a}{\sigma}\right) -
#'                                (x-\mu)\Phi\left(\frac{x-\mu-a}{\sigma}\right) +
#'                                \phi\left(\frac{x-\mu+a}{\sigma}\right) -
#'                                \phi\left(\frac{x-\mu-a}{\sigma}\right)\right]
#' }{
#' F(z) = \sigma/(2*a) * ((x-\mu)*\Phi((x-\mu+a)/\sigma) - (x-\mu)*\Phi((x-\mu-a)/\sigma) +
#'                        \phi((x-\mu+a)/\sigma) - \phi((x-\mu-a)/\sigma))
#' }
#'
#' @references
#' Bhattacharjee, G.P., Pandit, S.N.N., and Mohan, R. (1963).
#' Dimensional chains involving rectangular and normal error-distributions.
#' Technometrics, 5, 404-406.
#' 
#' @examples 
#' 
#' x <- rbhatt(1e5, 5, 3, 5)
#' hist(x, 100, freq = FALSE)
#' curve(dbhatt(x, 5, 3, 5), -20, 20, col = "red", add = TRUE)
#' hist(pbhatt(x, 5, 3, 5))
#' plot(ecdf(x))
#' curve(pbhatt(x, 5, 3, 5), -20, 20, col = "red", lwd = 2, add = TRUE)
#'
#' @name Bhattacharjee
#' @aliases Bhattacharjee
#' @aliases dbhatt
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dbhatt <- function(x, mu = 0, sigma = 1, a = sigma, log = FALSE) {
  cpp_dbhatt(x, mu, sigma, a, log[1L])
}


#' @rdname Bhattacharjee
#' @export

pbhatt <- function(q, mu = 0, sigma = 1, a = sigma, lower.tail = TRUE, log.p = FALSE) {
  cpp_pbhatt(q, mu, sigma, a, lower.tail[1L], log.p[1L])
}


#' @rdname Bhattacharjee
#' @export

rbhatt <- function(n, mu = 0, sigma = 1, a = sigma) {
  if (length(n) > 1) n <- length(n)
  cpp_rbhatt(n, mu, sigma, a)
}

