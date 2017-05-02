

#' Discrete Laplace distribution
#' 
#' Probability mass, distribution function and random generation
#' for the discrete Laplace distribution parametrized by location and scale.
#' 
#' @param x,q             vector of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param scale           scale parameter; \code{0 < scale < 1}.
#' @param location        location parameter.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#' 
#' @details 
#' 
#' If \eqn{U \sim \mathrm{Geometric}(1-p)}{U ~ Geometric(1-p)} and
#' \eqn{V \sim \mathrm{Geometric}(1-p)}{V ~ Geometric(1-p)},
#' then \eqn{U-V \sim \mathrm{DiscreteLaplace}(p)}{U-V ~ DiscreteLaplace(p)},
#' where geometric distribution is related to discrete Laplace distribution
#' in similar way as exponential distribution is related to Laplace distribution.
#' 
#' Probability mass function
#' 
#' \deqn{
#' f(x) = \frac{1-p}{1+p} p^{|x-\mu|}
#' }{
#' f(x) = (1-p)/(1+p) * p^(|x-\mu|)
#' }
#' 
#' Cumulative distribution function
#' 
#' \deqn{
#' F(x) = \left\{\begin{array}{ll}
#' \frac{p^{-|x-\mu|}}{1+p} & x < 0 \\
#' 1 - \frac{p^{|x-\mu|+1}}{1+p} & x \ge 0
#' \end{array}\right.
#' }{
#' F(x) = [if x < 0:] p^-floor(x-\mu))/(1+p) [else:] 1-(p^(floor(x-\mu)+1))/(1+p)
#' }
#' 
#' @references 
#' Inusah, S., & Kozubowski, T.J. (2006). A discrete analogue of the Laplace distribution.
#' Journal of statistical planning and inference, 136(3), 1090-1102.
#' 
#' @references 
#' Kotz, S., Kozubowski, T., & Podgorski, K. (2012).
#' The Laplace distribution and generalizations: a revisit with applications
#' to communications, economics, engineering, and finance.
#' Springer Science & Business Media.
#' 
#' @examples 
#' 
#' p <- 0.45
#' x <- rdlaplace(1e5, 0, p)
#' xx <- seq(-200, 200, by = 1)
#' plot(prop.table(table(x)))
#' lines(xx, ddlaplace(xx, 0, p), col = "red")
#' hist(pdlaplace(x, 0, p))
#' plot(ecdf(x))
#' lines(xx, pdlaplace(xx, 0, p), col = "red", type = "s")
#' 
#' @name DiscreteLaplace
#' @aliases DiscreteLaplace
#' @aliases ddlaplace
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Discrete
#' 
#' @export

ddlaplace <- function(x, location, scale, log = FALSE) {
  cpp_ddlaplace(x, location, scale, log[1L])
}


#' @rdname DiscreteLaplace
#' @export

pdlaplace <- function(q, location, scale, lower.tail = TRUE, log.p = FALSE) {
  cpp_pdlaplace(q, location, scale, lower.tail[1L], log.p[1L])
}


#' @rdname DiscreteLaplace
#' @export

rdlaplace <- function(n, location, scale) {
  if (length(n) > 1) n <- length(n)
  cpp_rdlaplace(n, location, scale)
}

