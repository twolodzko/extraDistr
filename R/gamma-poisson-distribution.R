

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
#' Probability density function (parametrized by scale)
#' \deqn{
#' f(x) = \frac{\Gamma(\alpha+x)}{x! \Gamma(\alpha)} \left(\frac{\beta}{1+\beta}\right)^k \left(1-\frac{\beta}{1+\beta}\right)^\alpha
#' }{
#' f(x) = \Gamma(\alpha+x) / (x!*\Gamma(\alpha)) * (\beta/(1+\beta))^x * (1-\beta/(1+\beta))^\alpha
#' }
#' 
#' \emph{Warning:} cumulative distribution function is defined as
#' \deqn{F(x) = \sum_{k=0}^x f(k)}{F(x) = f(0)+...+f(x)}
#' so it may be slow for large datasets.
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
#' plot(ecdf(x))
#' lines(xx, pgpois(xx, 7, 0.002), col = "red", lwd = 2)
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

pgpois <- function(q, shape, rate, scale = 1/rate, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_pgpois', PACKAGE = 'extraDistr', q, shape, scale, lower.tail, log.p)
}


#' @rdname GammaPoiss
#' @export

rgpois <- function(n, shape, rate, scale = 1/rate) {
  .Call('extraDistr_cpp_rgpois', PACKAGE = 'extraDistr', n, shape, scale)
}

