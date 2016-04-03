

#' Logarythmic Series distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Logarythmic Series distribution.
#'
#' @param x               matrix of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param theta           vector; concentration parameter.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability mass function
#' \deqn{
#' f(x) = \frac{-1}{\log(1-\theta)} \frac{\theta^x}{x}
#' }{
#' f(x) = (-1/log(1-\theta)*\theta^x) / x
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = \frac{-1}{\log(1-\theta)} \sum_{k=1}^x \frac{\theta^x}{x}
#' }{
#' F(x) = -1/log(1-\theta) * sum((\theta^x)/x)
#' }
#'
#' Quantile function and random generation are computed using
#' algorithm described in Krishnamoorthy (2006).
#'
#' @references
#' Krishnamoorthy, K. (2006). Handbook of Statistical Distributions
#' with Applications. Chapman & Hall/CRC
#'
#' @references
#' Forbes, C., Evans, M. Hastings, N., & Peacock, B. (2011).
#' Statistical Distributions. John Wiley & Sons.
#' 
#' @examples 
#' 
#' x <- rlgser(1e5, 0.66)
#' xx <- seq(0, 100, by = 1)
#' plot(prop.table(table(x)), type = "h")
#' lines(xx, dlgser(xx, 0.66), col = "red")
#' 
#' # Notice: distribution of F(X) is far from uniform:
#' hist(plgser(x, 0.66), 50)
#' 
#' plot(ecdf(x))
#' lines(xx, plgser(xx, 0.66), col = "red", lwd = 2)
#'
#' @name LogSeries
#' @aliases LogSeries
#' @aliases dlgser
#' @export

dlgser <- function(x, theta, log = FALSE) {
  .Call('extraDistr_cpp_dlgser', PACKAGE = 'extraDistr', x, theta, log)
}


#' @rdname LogSeries
#' @export

plgser <- function(x, theta, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_plgser', PACKAGE = 'extraDistr', x, theta, lower.tail, log.p)
}


#' @rdname LogSeries
#' @export

qlgser <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qlgser', PACKAGE = 'extraDistr', p, theta, lower.tail, log.p)
}


#' @rdname LogSeries
#' @export

rlgser <- function (n, theta) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rlgser', PACKAGE = 'extraDistr', n, theta)
}

