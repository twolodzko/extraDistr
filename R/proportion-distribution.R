

#' Beta distribution of proportions
#' 
#' Probability mass function, distribution function and random generation
#' for the reparametrized beta distribution.
#' 
#' @param x,q	            vector of quantiles.
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
#' Beta can be understood as a distribution of \eqn{x = k/n} proportions in
#' \eqn{n} trials where the average proportion is denoted as \eqn{\mu},
#' so it's parameters become \eqn{\alpha = n\mu+1} and
#' \eqn{\beta = n(1-\mu)+1} and it's density function becomes:
#' 
#' \deqn{
#' f(x) = \frac{1}{\mathrm{B}(n\mu+1, n(1-\mu)+1)} x^{n\mu} (1-x)^{n(1-\mu)}
#' }{
#' f(x) = 1/(B(n\mu+1, n(1-\mu)+1)) * x^(n\mu) * (1-x)^(n(1-\mu))
#' }
#' 
#' Alternatively \eqn{n} may be understood as precision parameter
#' as described Ferrari and Cribari-Neto (2004).
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

pprop <- function(q, size, mean, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_pprop', PACKAGE = 'extraDistr', q, size, mean, lower.tail, log.p)
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

