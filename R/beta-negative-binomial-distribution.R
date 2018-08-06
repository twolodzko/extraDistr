

#' Beta-negative binomial distribution
#'
#' Probability mass function and random generation
#' for the beta-negative binomial distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param alpha,beta      non-negative parameters of the beta distribution.
#' @param size            number of trials (zero or more). Must be strictly positive, need not be integer.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' If \eqn{p \sim \mathrm{Beta}(\alpha, \beta)}{p ~ Beta(\alpha, \beta)} and
#' \eqn{X \sim \mathrm{NegBinomial}(r, p)}{X ~ NegBinomial(r, p)}, then 
#' \eqn{X \sim \mathrm{BetaNegBinomial}(r, \alpha, \beta)}{X ~ BetaNegBinomial(r, \alpha, \beta)}.
#'
#' Probability mass function
#' \deqn{
#' f(x) = \frac{\Gamma(r+x)}{x! \,\Gamma(r)}
#'        \frac{\mathrm{B}(\alpha+r, \beta+x)}{\mathrm{B}(\alpha, \beta)}
#' }{
#' f(x) = \Gamma(r+x)/(x! \Gamma(r)) * B(\alpha+r, \beta+x) / B(\alpha, \beta)
#' }
#'
#' Cumulative distribution function is calculated using recursive algorithm that employs the fact that
#' \eqn{\Gamma(x) = (x - 1)!} and
#' \eqn{
#' \mathrm{B}(x, y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)}
#' }{
#' B(x, y) = (\Gamma(x)\Gamma(y))/\Gamma(x+y)
#' }. This enables re-writing probability mass function as
#' 
#' \deqn{
#' f(x) = \frac{(r+x-1)!}{x! \, \Gamma(r)} \frac{\frac{(\alpha+r-1)!\,(\beta+x-1)!}{(\alpha+\beta+r+x-1)!}}{\mathrm{B}(\alpha,\beta)}
#' }{
#' f(x) = ((r+x-1)!)/(x!*\Gamma(r))*(((\alpha+r-1)!*(\beta+x-1)!)/((\alpha+\beta+r+x-1)!))/B(\alpha,\beta)
#' }
#' 
#' what makes recursive updating from \eqn{x} to \eqn{x+1} easy using the properties of factorials
#' 
#' \deqn{
#' f(x+1) = \frac{(r+x-1)!\,(r+x)}{x!\,(x+1) \, \Gamma(r)} \frac{\frac{(\alpha+r-1)!\,(\beta+x-1)!\,(\beta+x)}{(\alpha+\beta+r+x-1)!\,(\alpha+\beta+r+x)}}{\mathrm{B}(\alpha,\beta)}
#' }{
#' f(x+1) = ((r+x-1)!*(r+x))/(x!*(x+1)*\Gamma(r))*(((\alpha+r-1)!*(\beta+x-1)!*(\beta+x))/((\alpha+\beta+r+x-1)!*(\alpha+\beta+r+x)))/B(\alpha,\beta)
#' }
#' 
#' and let's us efficiently calculate cumulative distribution function as a sum of probability mass functions
#' 
#' \deqn{F(x) = \sum_{k=0}^x f(k)}{F(x) = f(0)+...+f(x)}
#' 
#'
#' @seealso \code{\link[stats]{Beta}}, \code{\link[stats]{NegBinomial}}
#' 
#' @examples 
#' 
#' x <- rbnbinom(1e5, 1000, 5, 13)
#' xx <- 0:1e5
#' hist(x, 100, freq = FALSE)
#' lines(xx-0.5, dbnbinom(xx, 1000, 5, 13), col = "red")
#' hist(pbnbinom(x, 1000, 5, 13))
#' xx <- seq(0, 1e5, by = 0.1)
#' plot(ecdf(x))
#' lines(xx, pbnbinom(xx, 1000, 5, 13), col = "red", lwd = 2)
#'
#' @name BetaNegBinom
#' @aliases BetaNegBinom
#' @aliases dbnbinom
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Discrete
#'
#' @export

dbnbinom <- function(x, size, alpha = 1, beta = 1, log = FALSE) {
  cpp_dbnbinom(x, size, alpha, beta, log[1L])
}


#' @rdname BetaNegBinom
#' @export

pbnbinom <- function(q, size, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_pbnbinom(q, size, alpha, beta, lower.tail[1L], log.p[1L])
}


#' @rdname BetaNegBinom
#' @export

rbnbinom <- function(n, size, alpha = 1, beta = 1) {
  if (length(n) > 1) n <- length(n)
  cpp_rbnbinom(n, size, alpha, beta)
}

