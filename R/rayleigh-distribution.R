

#' Rayleigh distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Rayleigh distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param sigma           positive valued parameter.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = \frac{x}{\sigma^2} \exp\left(-\frac{x^2}{2\sigma^2}\right)
#' }{
#' f(x) = x/\sigma^2 * exp(-(x^2 / 2*\sigma^2))
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = 1 - \exp\left(-\frac{x^2}{2\sigma^2}\right)
#' }{
#' F(x) = 1 - exp(-x^2 / 2*\sigma^2)
#' }
#'
#' Quantile function
#' \deqn{
#' F^{-1}(p) = \sqrt{-2\sigma^2 \log(1-p)}
#' }{
#' F^-1(p) = sqrt(-2*\sigma^2 * log(1-p))
#' }
#'
#' @references
#' Krishnamoorthy, K. (2006). Handbook of Statistical Distributions
#' with Applications. Chapman & Hall/CRC.
#'
#' @references
#' Forbes, C., Evans, M. Hastings, N., & Peacock, B. (2011).
#' Statistical Distributions. John Wiley & Sons.
#' 
#' @examples 
#' 
#' x <- rrayleigh(1e5, 13)
#' hist(x, 100, freq = FALSE)
#' curve(drayleigh(x, 13), 0, 60, col = "red", add = TRUE)
#' hist(prayleigh(x, 13)) 
#' plot(ecdf(x))
#' curve(prayleigh(x, 13), 0, 60, col = "red", lwd = 2, add = TRUE) 
#'
#' @name Rayleigh
#' @aliases Rayleigh
#' @aliases drayleigh
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

drayleigh <- function(x, sigma = 1, log = FALSE) {
  cpp_drayleigh(x, sigma, log[1L])
}


#' @rdname Rayleigh
#' @export

prayleigh <- function(q, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_prayleigh(q, sigma, lower.tail[1L], log.p[1L])
}


#' @rdname Rayleigh
#' @export

qrayleigh <- function(p, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_qrayleigh(p, sigma, lower.tail[1L], log.p[1L])
}


#' @rdname Rayleigh
#' @export

rrayleigh <- function(n, sigma = 1) {
  if (length(n) > 1) n <- length(n)
  cpp_rrayleigh(n, sigma)
}
