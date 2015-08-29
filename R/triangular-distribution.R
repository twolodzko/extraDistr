

#' Triangular distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Triangular distribution.
#'
#' @param x 	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param a,b,c           minimum, maximum and shape parameters.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = \left\{\begin{array}{ll}
#' \frac{2(x-a)}{(b-a)(c-a)} & x < c \\
#' \frac{2}{b-a}             & x = c \\
#' \frac{2(b-x)}{(b-a)(b-c)} & x > c
#' \end{array}\right.
#' }{
#' f(x) = [if x < c:] (2*(x-a)) / ((b-a)*(c-a))
#'        [if x = c:] 2/(b-a)
#'        [if x >= c:] (2*(b-x)) / ((b-a)*(b-c))
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = \left\{\begin{array}{ll}
#' \frac{(x-a)^2}{(b-a)(c-a)}     & x \leq c \\
#' 1 - \frac{(b-x)^2}{(b-a)(b-c)} & x > c
#' \end{array}\right.
#' }{
#' F(x) = [if x <= c:] (x-a)^2 / ((b-a)*(c-a))
#'        [if x > c:] 1 - ((b-x)^2 / ((b-a)*(b-c)))
#' }
#'
#' Quantile function
#' \deqn{
#' F^{-1}(p) = \left\{\begin{array}{ll}
#' a + \sqrt{p \times (b-a)(c-a)} & p \leq \frac{c-a}{b-a} \\
#' b - \sqrt{(1-p)(b-a)(b-c)}     & p > \frac{c-a}{b-a}
#' \end{array}\right.
#' }{
#' F^-1(p) = [if p < (c-a)/(b-a):] a + sqrt(p*(b-a)*(c-a))
#'           [else:] b - sqrt((1-p)*(b-a)*(b-c))
#' }
#'
#' @references
#' Forbes, C., Evans, M. Hastings, N., & Peacock, B. (2011).
#' Statistical Distributions. John Wiley & Sons.
#'
#' @name Triangular
#' @aliases Triangular
#' @aliases dtriang
#' @keywords distribution
#'
#' @export

dtriang <- function(x, a, b, c = (a+b)/2, log = FALSE) {
  .Call('extraDistr_cpp_dtriang', PACKAGE = 'extraDistr', x, a, b, c, log)
}


#' @rdname Triangular
#' @export

ptriang <- function(x, a, b, c = (a+b)/2, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_ptriang', PACKAGE = 'extraDistr', x, a, b, c, lower.tail, log.p)
}


#' @rdname Triangular
#' @export

qtriang <- function(p, a, b, c = (a+b)/2, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qtriang', PACKAGE = 'extraDistr', p, a, b, c, lower.tail, log.p)
}


#' @rdname Triangular
#' @export

rtriang <- function(n, a, b, c = (a+b)/2) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rtriang', PACKAGE = 'extraDistr', n, a, b, c)
}

