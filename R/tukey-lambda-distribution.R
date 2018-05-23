

#' Tukey lambda distribution
#'
#' Quantile function, and random generation for the Tukey lambda
#' distribution.
#'
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param lambda	        shape parameter.   
#' @param log.p	          logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#' 
#' Tukey lambda distribution is a continuous probability distribution defined in terms
#' of its quantile function. It is typically used to identify other distributions.
#' 
#' Quantile function:
#' 
#' \deqn{F^{-1}(p) = \left\{\begin{array}{ll}
#' \frac{1}{\lambda} [p^\lambda - (1-p)^\lambda] & \lambda \ne 0 \\
#' \log(\frac{p}{1-p}) & \lambda = 0
#' \end{array}\right.
#' }{
#' F^-1(p) = [if \lambda != 0:] (p^\lambda - (1-p)^\lambda)/\lambda
#' [if \lambda = 0:] log(p/(1-p))
#' }
#'
#' @references
#' 
#' Joiner, B.L., & Rosenblatt, J.R. (1971).
#' Some properties of the range in samples from Tukey's symmetric lambda distributions.
#' Journal of the American Statistical Association, 66(334), 394-399.
#' 
#' @references 
#' 
#' Hastings Jr, C., Mosteller, F., Tukey, J.W., & Winsor, C.P. (1947).
#' Low moments for small samples: a comparative study of order statistics.
#' The Annals of Mathematical Statistics, 413-426.
#' 
#' @examples 
#' 
#' pp = seq(0, 1, by = 0.001)
#' partmp <- par(mfrow = c(2,3))
#' plot(qtlambda(pp, -1), pp, type = "l", main = "lambda = -1 (Cauchy)")
#' plot(qtlambda(pp, 0), pp, type = "l", main = "lambda = 0 (logistic)")
#' plot(qtlambda(pp, 0.14), pp, type = "l", main = "lambda = 0.14 (normal)")
#' plot(qtlambda(pp, 0.5), pp, type = "l", main = "lambda = 0.5 (concave)")
#' plot(qtlambda(pp, 1), pp, type = "l", main = "lambda = 1 (uniform)")
#' plot(qtlambda(pp, 2), pp, type = "l", main = "lambda = 2 (uniform)")
#' 
#' hist(rtlambda(1e5, -1), freq = FALSE, main = "lambda = -1 (Cauchy)")
#' hist(rtlambda(1e5, 0), freq = FALSE, main = "lambda = 0 (logistic)")
#' hist(rtlambda(1e5, 0.14), freq = FALSE, main = "lambda = 0.14 (normal)")
#' hist(rtlambda(1e5, 0.5), freq = FALSE, main = "lambda = 0.5 (concave)")
#' hist(rtlambda(1e5, 1), freq = FALSE, main = "lambda = 1 (uniform)")
#' hist(rtlambda(1e5, 2), freq = FALSE, main = "lambda = 2 (uniform)")
#' par(partmp)
#'
#' @name TukeyLambda
#' @aliases TukeyLambda
#' @aliases qtlambda
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

qtlambda <- function(p, lambda, lower.tail = TRUE, log.p = FALSE) {
  cpp_qtlambda(p, lambda, lower.tail[1L], log.p[1L])
}


#' @rdname TukeyLambda
#' @export

rtlambda <- function(n, lambda) {
  if (length(n) > 1) n <- length(n)
  cpp_rtlambda(n, lambda)
}

