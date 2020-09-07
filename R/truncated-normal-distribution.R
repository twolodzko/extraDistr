

#' Truncated normal distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the truncated normal distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param mean,sd         location and scale parameters. Scale must be positive.
#' @param a,b             lower and upper truncation points (\code{a < x <= b},
#'                        with \code{a = -Inf} and \code{b = Inf} by default).
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = \frac{\phi(\frac{x-\mu}{\sigma})}
#'             {\Phi(\frac{b-\mu}{\sigma}) - \Phi(\frac{a-\mu}{\sigma})}
#' }{
#' f(x) = \phi((x-\mu)/\sigma) / (\Phi((b-\mu)/\sigma) - \Phi((a-\mu)/\sigma))
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = \frac{\Phi(\frac{x-\mu}{\sigma}) - \Phi(\frac{a-\mu}{\sigma})}
#'             {\Phi(\frac{b-\mu}{\sigma}) - \Phi(\frac{a-\mu}{\sigma})}
#' }{
#' F(x) = (\Phi((x-\mu)/\sigma) - \Phi((a-\mu)/\sigma)) / (\Phi((b-\mu)/\sigma) - \Phi((a-\mu)/\sigma))
#' }
#'
#' Quantile function
#' \deqn{
#' F^{-1}(p) = \Phi^{-1}\left(\Phi\left(\frac{a-\mu}{\sigma}\right) + p \times
#'                      \left[\Phi\left(\frac{b-\mu}{\sigma}\right) -
#'                      \Phi\left(\frac{a-\mu}{\sigma}\right)\right]\right)
#' }{
#' F^-1(p) = \Phi^-1(\Phi((a-\mu)/\sigma) + p * (\Phi((b-\mu)/\sigma) - \Phi((a-\mu)/\sigma)))
#' }
#'
#' For random generation algorithm described by Robert (1995) is used.
#'
#' @references
#' Robert, C.P. (1995). Simulation of truncated normal variables.
#' Statistics and Computing 5(2): 121-125. \url{https://arxiv.org/abs/0907.4010}
#'
#' @references
#' Burkardt, J. (17 October 2014). The Truncated Normal Distribution. Florida State University.
#' \url{https://people.sc.fsu.edu/~jburkardt/presentations/truncated_normal.pdf}
#' 
#' @examples 
#' 
#' x <- rtnorm(1e5, 5, 3, b = 7)
#' hist(x, 100, freq = FALSE)
#' curve(dtnorm(x, 5, 3, b = 7), -8, 8, col = "red", add = TRUE)
#' hist(ptnorm(x, 5, 3, b = 7))
#' plot(ecdf(x))
#' curve(ptnorm(x, 5, 3, b = 7), -8, 8, col = "red", lwd = 2, add = TRUE)
#' 
#' R <- 1e5
#' partmp <- par(mfrow = c(2,4), mar = c(2,2,2,2))
#' 
#' hist(rtnorm(R), freq= FALSE, main = "", xlab = "", ylab = "")
#' curve(dtnorm(x), -5, 5, col = "red", add = TRUE)
#' 
#' hist(rtnorm(R, a = 0), freq= FALSE, main = "", xlab = "", ylab = "")
#' curve(dtnorm(x, a = 0), -1, 5, col = "red", add = TRUE)
#' 
#' hist(rtnorm(R, b = 0), freq= FALSE, main = "", xlab = "", ylab = "")
#' curve(dtnorm(x, b = 0), -5, 5, col = "red", add = TRUE)
#' 
#' hist(rtnorm(R, a = 0, b = 1), freq= FALSE, main = "", xlab = "", ylab = "")
#' curve(dtnorm(x, a = 0, b = 1), -1, 2, col = "red", add = TRUE)
#' 
#' hist(rtnorm(R, a = -1, b = 0), freq= FALSE, main = "", xlab = "", ylab = "")
#' curve(dtnorm(x, a = -1, b = 0), -2, 2, col = "red", add = TRUE)
#' 
#' hist(rtnorm(R, mean = -6, a = 0), freq= FALSE, main = "", xlab = "", ylab = "")
#' curve(dtnorm(x, mean = -6, a = 0), -2, 1, col = "red", add = TRUE)
#' 
#' hist(rtnorm(R, mean = 8, b = 0), freq= FALSE, main = "", xlab = "", ylab = "")
#' curve(dtnorm(x, mean = 8, b = 0), -2, 1, col = "red", add = TRUE)
#' 
#' hist(rtnorm(R, a = 3, b = 5), freq= FALSE, main = "", xlab = "", ylab = "")
#' curve(dtnorm(x, a = 3, b = 5), 2, 5, col = "red", add = TRUE)
#' 
#' par(partmp)
#'
#' @name TruncNormal
#' @aliases TruncNormal
#' @aliases dtnorm
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dtnorm <- function(x, mean = 0, sd = 1, a = -Inf, b = Inf, log = FALSE) {
  cpp_dtnorm(x, mean, sd, a, b, log[1L])
}


#' @rdname TruncNormal
#' @export

ptnorm <- function(q, mean = 0, sd = 1, a = -Inf, b = Inf, lower.tail = TRUE, log.p = FALSE) {
  cpp_ptnorm(q, mean, sd, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname TruncNormal
#' @export

qtnorm <- function(p, mean = 0, sd = 1, a = -Inf, b = Inf, lower.tail = TRUE, log.p = FALSE) {
  cpp_qtnorm(p, mean, sd, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname TruncNormal
#' @export

rtnorm <- function(n, mean = 0, sd = 1, a = -Inf, b = Inf) {
  if (length(n) > 1) n <- length(n)
  cpp_rtnorm(n, mean, sd, a, b)
}

