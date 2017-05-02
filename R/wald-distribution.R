

#' Wald (inverse Gaussian) distribution
#'
#' Density, distribution function and random generation
#' for the Wald distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param mu,lambda	      location and shape parameters. Scale must be positive.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = \sqrt{\frac{\lambda}{2\pi x^3}} \exp\left( \frac{-\lambda(x-\mu)^2}{2\mu^2 x} \right)
#' }{
#' f(x) = sqrt(\lambda/(2*\pi*x^3)) * exp((-\lambda*(x-\mu)^2)/(2*\mu^2*x))
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = \Phi\left(\sqrt{\frac{\lambda}{x}} \left(\frac{x}{\mu}-1 \right) \right) +
#' \exp\left(\frac{2\lambda}{\mu} \right) \Phi\left(\sqrt{\frac{\lambda}{x}}
#' \left(\frac{x}{\mu}+1 \right) \right)
#' }{
#' F(x) = \Phi(sqrt(\lambda/\mu)*(x/\mu-1)) - exp((2*\lambda)/\mu) *
#' \Phi(sqrt(\lambda/\mu)*(x/\mu+1))
#' }
#' 
#' Random generation is done using the algorithm described by Michael, Schucany and Haas (1976).
#' 
#' @references 
#' 
#' Michael, J.R., Schucany, W.R., and Haas, R.W. (1976).
#' Generating Random Variates Using Transformations with Multiple Roots.
#' The American Statistician, 30(2): 88-90.
#' 
#' @examples 
#' 
#' x <- rwald(1e5, 5, 16)
#' hist(x, 100, freq = FALSE)
#' curve(dwald(x, 5, 16), 0, 50, col = "red", add = TRUE)
#' hist(pwald(x, 5, 16))
#' plot(ecdf(x))
#' curve(pwald(x, 5, 16), 0, 50, col = "red", lwd = 2, add = TRUE)
#'
#' @name Wald
#' @aliases Wald
#' @aliases dwald
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dwald <- function(x, mu, lambda, log = FALSE) {
  cpp_dwald(x, mu, lambda, log[1L])
}


#' @rdname Wald
#' @export

pwald <- function(q, mu, lambda, lower.tail = TRUE, log.p = FALSE) {
  cpp_pwald(q, mu, lambda, lower.tail[1L], log.p[1L])
}


#' @rdname Wald
#' @export

rwald <- function(n, mu, lambda) {
  if (length(n) > 1) n <- length(n)
  cpp_rwald(n, mu, lambda)
}

