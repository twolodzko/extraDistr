

#' Discrete gamma distribution
#' 
#' Probability mass function, distribution function and random generation
#' for discrete gamma distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param rate	          an alternative way to specify the scale.
#' @param shape,scale	    shape and scale parameters. Must be positive, scale strictly.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'                
#' @details 
#'                                 
#' Probability mass function of discrete gamma distribution \eqn{f_Y(y)}{f}
#' is defined by discretization of continuous gamma distribution
#' \eqn{f_Y(y) = S_X(y) - S_X(y+1)}{f(y) = S(x) - S(x+1)}
#' where \eqn{S_X}{S} is a survival function of continuous gamma distribution.
#' 
#' @references 
#' Chakraborty, S. and Chakravarty, D. (2012).
#' Discrete Gamma distributions: Properties and parameter estimations.
#' Communications in Statistics-Theory and Methods, 41(18), 3301-3324.
#' 
#' @seealso \code{\link[stats]{GammaDist}}, \code{\link{DiscreteNormal}}
#' 
#' @examples 
#' 
#' x <- rdgamma(1e5, 9, 1)
#' xx <- 0:50
#' plot(prop.table(table(x)))
#' lines(xx, ddgamma(xx, 9, 1), col = "red")
#' hist(pdgamma(x, 9, 1))
#' plot(ecdf(x))
#' xx <- seq(0, 50, 0.1)
#' lines(xx, pdgamma(xx, 9, 1), col = "red", lwd = 2, type = "s")
#' 
#' @name DiscreteGamma
#' @aliases DiscreteGamma
#' @aliases ddgamma
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Discrete
#' 
#' @export

ddgamma <- function(x, shape, rate = 1, scale = 1/rate, log = FALSE) {
  if (!missing(rate) && !missing(scale)) {
    if (abs(rate * scale - 1) < 1e-15)
      warning("specify 'rate' or 'scale' but not both")
    else stop("specify 'rate' or 'scale' but not both")
  }
  cpp_ddgamma(x, shape, scale, log[1L])
}


#' @rdname DiscreteGamma
#' @export

pdgamma <- function(q, shape, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE) {
  pgamma(floor(q)+1, shape, scale = scale, lower.tail = lower.tail[1L], log.p = log.p[1L])
}


#' @rdname DiscreteGamma
#' @export

rdgamma <- function(n, shape, rate = 1, scale = 1/rate) {
  floor(rgamma(n, shape, scale = scale))
}

