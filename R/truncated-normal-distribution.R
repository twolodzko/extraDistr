

#' Truncated Normal distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Truncated Normal distribution.
#'
#' @param x 	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param mu,sigma        location and scale parameters.
#' @param a,b             minimal and maximal boundries for truncation
#'                        (\code{-Inf} and \code{Inf} by default).
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = \frac{\phi(z)}
#'             {\Phi(\frac{b-\mu}{\sigma}) - \Phi(\frac{a-\mu}{\sigma})}
#' }{
#' f(x) = \phi(z) / (\Phi((b-\mu)/\sigma) - \Phi((a-\mu)/\sigma))
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = \frac{\Phi(z) - \Phi(\frac{a-\mu}{\sigma})}
#'             {\Phi(\frac{b-\mu}{\sigma}) - \Phi(\frac{a-\mu}{\sigma})}
#' }{
#' F(x) = (\Phi(z) - \Phi((a-\mu)/\sigma)) / (\Phi((b-\mu)/\sigma) - \Phi((a-\mu)/\sigma))
#' }
#'
#' Quantile function
#' \deqn{
#' F^{-1}(p) = \Phi^{-1}\left(\Phi(\frac{a-\mu}{\sigma}) + p \times
#'                      (\Phi(\frac{b-\mu}{\sigma}) - \Phi(\frac{a-\mu}{\sigma}))\right)
#' }{
#' F^-1(p) = \Phi^-1(\Phi((a-\mu)/\sigma) + p * (\Phi((b-\mu)/\sigma) - \Phi((a-\mu)/\sigma)))
#' }
#'
#' where \eqn{ z = \frac{x-\mu}{\sigma} }{ z = (x-\mu)/\sigma }.
#'
#' For random generation algorithm described by Robert (1995) is used.
#'
#' @references
#' Robert, C.P. (1995). Simulation of truncated normal variables.
#' Statistics and Computing 5(2): 121-125. \url{http://arxiv.org/abs/0907.4010}
#'
#' @references
#' Burkardt, J. (17 October 2014). The Truncated Normal Distribution. Florida State University.
#' \url{http://people.sc.fsu.edu/~jburkardt/presentations/truncated_normal.pdf}
#'
#' @name TruncNormal
#' @aliases TruncNormal
#' @aliases dtnorm
#' @keywords distribution
#'
#' @export

dtnorm <- function(x, mu = 0, sigma = 1, a = -Inf, b = Inf, log = FALSE) {
  .Call('extraDistr_cpp_dtnorm', PACKAGE = 'extraDistr', x, mu, sigma, a, b, log)
}


#' @rdname TruncNormal
#' @export

ptnorm <- function(x, mu = 0, sigma = 1, a = -Inf, b = Inf, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_ptnorm', PACKAGE = 'extraDistr', x, mu, sigma, a, b, lower.tail, log.p)
}


#' @rdname TruncNormal
#' @export

qtnorm <- function(p, mu = 0, sigma = 1, a = -Inf, b = Inf, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qtnorm', PACKAGE = 'extraDistr', p, mu, sigma, a, b, lower.tail, log.p)
}


#' @rdname TruncNormal
#' @export

rtnorm <- function(n, mu = 0, sigma = 1, a = -Inf, b = Inf) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rtnorm', PACKAGE = 'extraDistr', n, mu, sigma, a, b)
}

