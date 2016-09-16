

#' Beta-negative binomial distribution
#'
#' Probability mass function and random generation
#' for the beta-negative binomial distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param alpha,beta      non-negative parameters of the beta distribution.
#' @param size            number of trials (zero or more).
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
#' f(x) = \frac{\Gamma(r+x)}{x! \Gamma(r)}
#'        \frac{\mathrm{B}(\alpha+r, \beta+x)}{\mathrm{B}(\alpha, \beta)}
#' }{
#' f(x) = \Gamma(r+x)/(x! \Gamma(r)) * B(\alpha+r, \beta+x) / B(\alpha, \beta)
#' }
#'
#' \emph{Warning:} cumulative distribution function is defined as
#' \deqn{F(x) = \sum_{k=0}^x f(k)}{F(x) = f(0)+...+f(x)}
#' so it may be slow for large datasets.
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
#' plot(ecdf(x))
#' lines(xx, pbnbinom(xx, 1000, 5, 13), col = "red", lwd = 2)
#'
#' @name BetaNegBinom
#' @aliases BetaNegBinom
#' @aliases dbnbinom
#' @keywords distribution
#'
#' @export

dbnbinom <- function(x, size, alpha = 1, beta = 1, log = FALSE) {
  cpp_dbnbinom(x, size, alpha, beta, log)
}


#' @rdname BetaNegBinom
#' @export

pbnbinom <- function(q, size, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_pbnbinom(q, size, alpha, beta, lower.tail, log.p)
}


#' @rdname BetaNegBinom
#' @export

rbnbinom <- function(n, size, alpha = 1, beta = 1) {
  if (length(n) > 1) n <- length(n)
  cpp_rbnbinom(n, size, alpha, beta)
}

