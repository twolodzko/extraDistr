

#' Half-Cauchy distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the half-Cauchy distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param sigma	          positive valued scale parameter.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'                        
#' @details
#' If \eqn{X} follows Cauchy centered at 0 and parametrized by scale \eqn{\sigma},
#' then \eqn{|X|} follows half-Cauchy distribution parametrized by
#' scale \eqn{\sigma}. Half-Cauchy distribution is a special case of half-t
#' distribution with \eqn{\nu=1} degrees of freedom.
#'                        
#' @references
#' Gelman, A. (2006). Prior distributions for variance parameters in hierarchical
#' models (comment on article by Browne and Draper).
#' Bayesian analysis, 1(3), 515-534.
#' 
#' @references 
#' Jacob, E. and Jayakumar, K. (2012).
#' On Half-Cauchy Distribution and Process.
#' International Journal of Statistika and Mathematika, 3(2), 77-81.
#' 
#' @seealso \code{\link{HalfT}}
#' 
#' @examples 
#' 
#' x <- rhcauchy(1e5, 2)
#' hist(x, 2e5, freq = FALSE, xlim = c(0, 100))
#' curve(dhcauchy(x, 2), 0, 100, col = "red", add = TRUE)
#' hist(phcauchy(x, 2))
#' plot(ecdf(x), xlim = c(0, 100))
#' curve(phcauchy(x, 2), col = "red", lwd = 2, add = TRUE)
#'
#' @name HalfCauchy
#' @aliases HalfCauchy
#' @aliases dhcauchy
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dhcauchy <- function(x, sigma = 1, log = FALSE) {
  cpp_dhcauchy(x, sigma, log[1L])
}


#' @rdname HalfCauchy
#' @export

phcauchy <- function(q, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_phcauchy(q, sigma, lower.tail[1L], log.p[1L])
}


#' @rdname HalfCauchy
#' @export

qhcauchy <- function(p, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_qhcauchy(p, sigma, lower.tail[1L], log.p[1L])
}


#' @rdname HalfCauchy
#' @export

rhcauchy <- function(n, sigma = 1) {
  if (length(n) > 1) n <- length(n)
  cpp_rhcauchy(n, sigma)
}

