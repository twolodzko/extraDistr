

#' "Huber density" distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the "Huber density" distribution.
#'
#' @param x,q	             vector of quantiles.
#' @param p	               vector of probabilities.
#' @param n	               number of observations. If \code{length(n) > 1},
#'                         the length is taken to be the number required.
#' @param mu,sigma,epsilon location, and scale, and shape parameters.
#'                         Scale and shape must be positive.
#' @param log,log.p	       logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	     logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                         otherwise, \eqn{P[X > x]}.
#'
#' @details
#' 
#' Huber density is connected to Huber loss and can be defined as:
#'
#' \deqn{
#' f(x) = \frac{1-\epsilon}{\sqrt{2\pi}} e^{-\rho_k(x)}
#' }{
#' f(x) = (1-\epsilon)/sqrt(2*\pi) * exp(-\rho(x, k))
#' }
#'
#' where
#'
#' \deqn{
#' \rho_k(x) =
#' \left\{\begin{array}{ll}
#' \frac{1}{2} x^2       & |x|\le k \\
#' k|x|- \frac{1}{2} k^2 & |x|>k
#' \end{array}\right.
#' }{
#' \rho(x, k) = [if abs(x) <= k:] (x^2)/2 [else:] k*abs(x) - (k^2)/2
#' }
#'
#' and \eqn{\epsilon} satisfies
#'
#' \deqn{
#' \frac{2\phi(k)}{k} - 2\Phi(-k) = \frac{\epsilon}{1-\epsilon}
#' }{
#' (2 * \phi(k))/k - 2* \Phi(-k) = \epsilon/(1-\epsilon)
#' }
#' 
#' @references
#' Huber, P.J. (1964). Robust Estimation of a Location Parameter.
#' Annals of Statistics, 53(1), 73-101.
#' 
#' @references
#' Huber, P.J. (1981). Robust Statistics. Wiley.
#' 
#' @references
#' Schumann, D. (2009). Robust Variable Selection. ProQuest.
#' 
#' @examples 
#' 
#' x <- rhuber(1e5, 5, 2, 3)
#' hist(x, 100, freq = FALSE)
#' curve(dhuber(x, 5, 2, 3), -20, 20, col = "red", add = TRUE)
#' hist(phuber(x, 5, 2, 3))
#' plot(ecdf(x))
#' curve(phuber(x, 5, 2, 3), -20, 20, col = "red", lwd = 2, add = TRUE)
#'
#' @name Huber
#' @aliases Huber
#' @aliases dhuber
#' @keywords distribution
#'
#' @export

dhuber <- function(x, mu = 0, sigma = 1, epsilon = 1.345, log = FALSE) {
  cpp_dhuber(x, mu, sigma, epsilon, log)
}


#' @rdname Huber
#' @export

phuber <- function(q, mu = 0, sigma = 1, epsilon = 1.345, lower.tail = TRUE, log.p = FALSE) {
  cpp_phuber(q, mu, sigma, epsilon, lower.tail, log.p)
}


#' @rdname Huber
#' @export

qhuber <- function(p, mu = 0, sigma = 1, epsilon = 1.345, lower.tail = TRUE, log.p = FALSE) {
  cpp_qhuber(p, mu, sigma, epsilon, lower.tail, log.p)
}


#' @rdname Huber
#' @export

rhuber <- function(n, mu = 0, sigma = 1, epsilon = 1.345) {
  if (length(n) > 1) n <- length(n)
  cpp_rhuber(n, mu, sigma, epsilon)
}

