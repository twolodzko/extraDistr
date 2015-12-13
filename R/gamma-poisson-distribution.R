

#' Gamma-Poisson distribution
#'
#' Probability mass function and random generation
#' for the Gamma-Poisson distribution.
#'
#' @param x 	            vector of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param rate	          an alternative way to specify the scale.
#' @param shape,scale	    shape and scale parameters. Must be positive,
#'                        scale strictly.
#' @param log     	      logical; if TRUE, probabilities p are given as log(p).
#'
#' @details
#' Gamma-Poisson distribution arises as a continuous mixture of
#' Poisson distributions, where where the mixing distribution
#' of the Poisson rate \eqn{\lambda} is a gamma distribution.
#' When \eqn{X \sim \mathrm{Poisson}(\lambda)}{X ~ Poisson(\lambda)}
#' and \eqn{\lambda \sim \mathrm{Gamma}(\alpha, \beta)}{\lambda ~ Gamma(\alpha, \beta)}, then \eqn{X \sim \mathrm{GammaPoisson}(\alpha, \beta)}{X ~ Gamma-Poisson(\alpha, \beta)}.
#'
#' Probability density function (parametrized by scale)
#' \deqn{
#' f(x) = \frac{\Gamma(\alpha+x)}{x! \Gamma(\alpha)} p^k (1-p)^\alpha
#' }{
#' f(x) = \Gamma(\alpha+x) / (x!*\Gamma(\alpha)) * p^x * (1-p)^\alpha
#' }
#' 
#' where \eqn{p = \frac{\beta}{1+\beta}}{p = \beta/(1+\beta)}.
#'
#' @seealso \code{\link{Gamma}}, \code{\link{Poisson}}
#'
#' @name GammaPoiss
#' @aliases GammaPoiss
#' @aliases dgpois
#' @keywords distribution
#'
#' @export

dgpois <- function(x, shape, rate, scale = 1/rate, log = FALSE) {
  .Call('extraDistr_cpp_dgpois', PACKAGE = 'extraDistr', x, shape, scale, log)
}


#' @rdname GammaPoiss
#' @export

rgpois <- function(n, shape, rate, scale = 1/rate) {
  .Call('extraDistr_cpp_rgpois', PACKAGE = 'extraDistr', n, shape, scale)
}

