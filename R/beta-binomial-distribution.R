

#' Beta-binomial distribution
#'
#' Probability mass function and random generation
#' for the beta-binomial distribution.
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
#' \eqn{X \sim \mathrm{Binomial}(n, p)}{X ~ Binomial(n, p)}, then 
#' \eqn{X \sim \mathrm{BetaBinomial}(n, \alpha, \beta)}{X ~ BetaBinomial(n, \alpha, \beta)}.
#' 
#' Probability mass function
#' \deqn{
#' f(x) = {n \choose x} \frac{\mathrm{B}(x+\alpha, n-x+\beta)}{\mathrm{B}(\alpha, \beta)}
#' }{
#' f(x) = choose(n, x) * B(x+\alpha, n-x+\beta) / B(\alpha, \beta)
#' }
#'
#' Cumulative distribution function is defined as
#' 
#' \deqn{F(x) = \sum_{k=0}^x f(k)}{F(x) = f(0)+...+f(x)}
#' 
#' and it is calculated using recursive algorithm that employs the fact that
#' \eqn{
#' \mathrm{B}(x, y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)}
#' }{
#' B(x, y) = (\Gamma(x)\Gamma(y))/\Gamma(x+y)
#' },
#' and \eqn{\Gamma(x) = (x - 1)!}, and
#' \eqn{
#' {n \choose k} = \prod_{i=1}^k \frac{n+1-i}{i}
#' }{
#' choose(n, k) = prod((n+1-(1:k))/(1:k))
#' }.
#'
#' @seealso \code{\link[stats]{Beta}}, \code{\link[stats]{Binomial}}
#' 
#' @examples 
#' 
#' x <- rbbinom(1e5, 1000, 5, 13)
#' xx <- 0:1000
#' hist(x, 100, freq = FALSE)
#' lines(xx-0.5, dbbinom(xx, 1000, 5, 13), col = "red")
#' hist(pbbinom(x, 1000, 5, 13))
#' plot(ecdf(x))
#' lines(xx, pbbinom(xx, 1000, 5, 13), col = "red", lwd = 2)
#'
#' @name BetaBinom
#' @aliases BetaBinom
#' @aliases dbbinom
#' @keywords distribution
#'
#' @export

dbbinom <- function(x, size, alpha = 1, beta = 1, log = FALSE) {
  cpp_dbbinom(x, size, alpha, beta, log)
}


#' @rdname BetaBinom
#' @export

pbbinom <- function(q, size, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_pbbinom(q, size, alpha, beta, lower.tail, log.p)
}


#' @rdname BetaBinom
#' @export

rbbinom <- function(n, size, alpha = 1, beta = 1) {
  if (length(n) > 1) n <- length(n)
  cpp_rbbinom(n, size, alpha, beta)
}

