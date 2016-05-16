

#' Beta distribution of proportions
#' 
#' Probability mass function, distribution function and random generation
#' for the reparametrized beta distribution.
#' 
#' @param x 	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param size            precision or number of binomial trials (zero or more).
#' @param mean            mean proportion or probability of success on each trial;
#'                        \code{0 < mean < 1}.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'                        
#' @details
#' 
#' Probability mass function of binomial distribution is
#' 
#' \deqn{
#' {n \choose k} p^k (1-p) ^{n-k} 
#' }{
#' choose(n,k) * p^k (1-p)^(n-k)
#' }
#' 
#' probability density function of beta distribution is
#' 
#' \deqn{
#' \frac{1}{\mathrm{B}(\alpha, \beta)} p^{\alpha-1} (1-p)^{\beta-1}
#' }{
#' 1/beta(\alpha, \beta) p^(\alpha-1) (1-p)^(\beta-1)
#' }
#' 
#' we can rewrite
#' 
#' \deqn{
#' {n \choose k} = \frac{1}{(n+1) \mathrm{B}(k+1, n-k+1)}
#' }{
#' choose(n,k) = 1/((n+1) * beta(k+1, n-k+1))
#' }
#' 
#' if we substitute \eqn{k+1 = \alpha} and \eqn{n-k+1 = \beta} then pmf
#' of binomial distribution becomes
#' 
#' \deqn{
#' \frac{1}{(n+1) \mathrm{B}(\alpha, \beta)} p^{\alpha-1} (1-p)^{\beta-1}
#' }{
#' 1/((n+1) * beta(\alpha, \beta)) * p^(\alpha-1) * (1-p)^(\beta-1)
#' }
#' 
#' so beta can be understood as a distribution of \eqn{k/n} proportions in
#' \eqn{n} trials where the average proportion is denoted as \eqn{\mu}
#' 
#' \deqn{
#' \frac{1}{\mathrm{B}(n\mu, n(1-\mu))} p^{n\mu+1} (1-p)^{n(1-\mu)+1}
#' }{
#' 1/(beta(n\mu, n(1-\mu))) * p^(n\mu+1) * (1-p)^(n(1-\mu)+1)
#' }
#' 
#' @references
#' Ferrari, S., & Cribari-Neto, F. (2004). Beta regression for modelling rates and proportions.
#' Journal of Applied Statistics, 31(7), 799-815.
#' 
#' @examples 
#' 
#' x <- rprop(1e5, 100, 0.33)
#' xx <- seq(0, 1, by = 0.01)
#' hist(x, 100, freq = FALSE)
#' lines(xx, dprop(xx, 100, 0.33), col = "red")
#' hist(pprop(x, 100, 0.33))
#' plot(ecdf(x))
#' lines(xx, pprop(xx, 100, 0.33), col = "red", lwd = 2)
#'                        
#' @name PropBeta
#' @aliases PropBeta
#' @aliases dprop
#' @keywords distribution
#'
#' @export

dprop <- function(x, size, mean, log = FALSE) {
  .Call('extraDistr_cpp_dprop', PACKAGE = 'extraDistr', x, size, mean, log)
}


#' @rdname PropBeta
#' @export

pprop <- function(x, size, mean, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_pprop', PACKAGE = 'extraDistr', x, size, mean, lower.tail, log.p)
}


#' @rdname PropBeta
#' @export

qprop <- function(p, size, mean, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qprop', PACKAGE = 'extraDistr', p, size, mean, lower.tail, log.p)
}


#' @rdname PropBeta
#' @export

rprop <- function(n, size, mean) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rprop', PACKAGE = 'extraDistr', n, size, mean)
}

