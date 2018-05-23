

#' Gamma-Poisson distribution
#'
#' Probability mass function and random generation
#' for the gamma-Poisson distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param rate	          an alternative way to specify the scale.
#' @param shape,scale	    shape and scale parameters. Must be positive,
#'                        scale strictly.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#' Gamma-Poisson distribution arises as a continuous mixture of
#' Poisson distributions, where the mixing distribution
#' of the Poisson rate \eqn{\lambda} is a gamma distribution.
#' When \eqn{X \sim \mathrm{Poisson}(\lambda)}{X ~ Poisson(\lambda)}
#' and \eqn{\lambda \sim \mathrm{Gamma}(\alpha, \beta)}{\lambda ~ Gamma(\alpha, \beta)}, then \eqn{X \sim \mathrm{GammaPoisson}(\alpha, \beta)}{X ~ Gamma-Poisson(\alpha, \beta)}.
#'
#' Probability mass function
#' \deqn{
#' f(x) = \frac{\Gamma(\alpha+x)}{x! \, \Gamma(\alpha)} \left(\frac{\beta}{1+\beta}\right)^x \left(1-\frac{\beta}{1+\beta}\right)^\alpha
#' }{
#' f(x) = \Gamma(\alpha+x) / (x!*\Gamma(\alpha)) * (\beta/(1+\beta))^x * (1-\beta/(1+\beta))^\alpha
#' }
#' 
#' Cumulative distribution function is calculated using recursive algorithm that employs the fact that
#' \eqn{\Gamma(x) = (x - 1)!}. This enables re-writing probability mass function as
#' 
#' \deqn{
#' f(x) = \frac{(\alpha+x-1)!}{x! \, \Gamma(\alpha)} \left( \frac{\beta}{1+\beta} \right)^x \left( 1- \frac{\beta}{1+\beta} \right)^\alpha
#' }{
#' f(x) = ((\alpha+x-1)!)/(x!*\Gamma(\alpha))*(\beta/(1+\beta))^x*(1-\beta/(1+\beta))^\alpha
#' }
#' 
#' what makes recursive updating from \eqn{x} to \eqn{x+1} easy using the properties of factorials
#' 
#' \deqn{
#' f(x+1) = \frac{(\alpha+x-1)! \, (\alpha+x)}{x! \,(x+1) \, \Gamma(\alpha)} \left( \frac{\beta}{1+\beta} \right)^x \left( \frac{\beta}{1+\beta} \right) \left( 1- \frac{\beta}{1+\beta} \right)^\alpha
#' }{
#' f(x+1) = ((\alpha+x-1)!*(\alpha+x))/(x!*(x+1)*\Gamma(\alpha))*(\beta/(1+\beta))^x*(\beta/(1+\beta))*(1-\beta/(1+\beta))^\alpha
#' }
#' 
#' and let's us efficiently calculate cumulative distribution function as a sum of probability mass functions
#' 
#' \deqn{F(x) = \sum_{k=0}^x f(k)}{F(x) = f(0)+...+f(x)}
#' 
#'
#' @seealso \code{\link[stats]{Gamma}}, \code{\link[stats]{Poisson}}
#' 
#' @examples 
#' 
#' x <- rgpois(1e5, 7, 0.002)
#' xx <- seq(0, 12000, by = 1)
#' hist(x, 100, freq = FALSE)
#' lines(xx, dgpois(xx, 7, 0.002), col = "red")
#' hist(pgpois(x, 7, 0.002))
#' xx <- seq(0, 12000, by = 0.1)
#' plot(ecdf(x))
#' lines(xx, pgpois(xx, 7, 0.002), col = "red", lwd = 2)
#'
#' @name GammaPoiss
#' @aliases GammaPoiss
#' @aliases dgpois
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Discrete
#'
#' @export

dgpois <- function(x, shape, rate, scale = 1/rate, log = FALSE) {
  cpp_dgpois(x, shape, scale, log[1L])
}


#' @rdname GammaPoiss
#' @export

pgpois <- function(q, shape, rate, scale = 1/rate, lower.tail = TRUE, log.p = FALSE) {
  cpp_pgpois(q, shape, scale, lower.tail[1L], log.p[1L])
}


#' @rdname GammaPoiss
#' @export

rgpois <- function(n, shape, rate, scale = 1/rate) {
  cpp_rgpois(n, shape, scale)
}

