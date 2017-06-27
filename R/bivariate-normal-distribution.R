

#' Bivariate normal distribution
#'
#' Density, distribution function and random generation
#' for the bivariate normal distribution.
#'
#' @param x,y	        vectors of quantiles; alternatively x may be a two-column
#'                    matrix (or data.frame) and y may be omitted.
#' @param n	          number of observations. If \code{length(n) > 1},
#'                    the length is taken to be the number required.
#' @param mean1,mean2 vectors of means.
#' @param sd1,sd2     vectors of standard deviations.
#' @param cor         vector of correlations (\code{-1 < cor < 1}).
#' @param log     	  logical; if TRUE, probabilities p are given as log(p).
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = \frac{1}{2\pi\sqrt{1-\rho^2}\sigma_1\sigma_2}
#'        \exp\left\{-\frac{1}{2(1-\rho^2)} \left[\left(\frac{x_1 - \mu_1}{\sigma_1}\right)^2 -
#'        2\rho \left(\frac{x_1 - \mu_1}{\sigma_1}\right) \left(\frac{x_2 - \mu_2}{\sigma_2}\right) +
#'        \left(\frac{x_2 - \mu_2}{\sigma_2}\right)^2\right]\right\}
#' }{
#' f(x) = 1/(2*\pi*sqrt(1-\rho^2)*\sigma1*\sigma2) * exp(-(1/(2*(1-\rho^2)*
#'        (((x1-\mu1)/\sigma1)^2 - 2*\rho*((x1-\mu1)/\sigma2)*((x2-\mu2)/\sigma2) *
#'        ((x2-\mu2)/\sigma2)^2))))
#' }
#'
#' @references
#' Krishnamoorthy, K. (2006). Handbook of Statistical Distributions
#' with Applications. Chapman & Hall/CRC
#' 
#' @references 
#' Mukhopadhyay, N. (2000). Probability and statistical inference.
#' Chapman & Hall/CRC
#' 
#' @examples 
#' 
#' y <- x <- seq(-4, 4, by = 0.25)
#' z <- outer(x, y, function(x, y) dbvnorm(x, y, cor = -0.75))
#' persp(x, y, z)
#'
#' y <- x <- seq(-4, 4, by = 0.25)
#' z <- outer(x, y, function(x, y) dbvnorm(x, y, cor = -0.25))
#' persp(x, y, z)
#'
#' @seealso \code{\link[stats]{Normal}}
#'
#' @name BivNormal
#' @aliases BivNormal
#' @aliases dbvnorm
#' 
#' @keywords distribution
#' @concept Multivariate
#' @concept Continuous
#'
#' @export

dbvnorm <- function(x, y = NULL, mean1 = 0, mean2 = mean1, sd1 = 1, sd2 = sd1, cor = 0, log = FALSE) {
  if (is.null(y)) {
    if ((is.matrix(x) || is.data.frame(x)) && ncol(x) == 2) {
      y <- x[, 2]
      x <- x[, 1]
    } else {
      stop("y is not provided while x is not a two-column matrix")
    }
  }
  cpp_dbnorm(x, y, mean1, mean2, sd1, sd2, cor, log[1L])
}


#' @rdname BivNormal
#' @export

rbvnorm <- function(n, mean1 = 0, mean2 = mean1, sd1 = 1, sd2 = sd1, cor = 0) {
  if (length(n) > 1) n <- length(n)
  cpp_rbnorm(n, mean1, mean2, sd1, sd2, cor)
}

