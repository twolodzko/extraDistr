

#' Mixture of Poisson distributions
#'
#' Density, distribution function and random generation
#' for the mixture of Poisson distributions.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param lambda 	        matrix (or vector) of (non-negative) means.
#' @param alpha           matrix (or vector) of mixing proportions;
#'                        mixing proportions need to sum up to 1.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = \alpha_1 f_1(x; \lambda_1) + \dots + \alpha_k f_k(x; \lambda_k)
#' }{
#' f(x) = \alpha[1] * f1(x; \lambda[1]) + \dots + \alpha[k] * fk(x; \lambda[k])
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = \alpha_1 F_1(x; \lambda_1) + \dots + \alpha_k F_k(x; \lambda_k)
#' }{
#' F(x) = \alpha[1] * F1(x; \lambda[1]) + \dots + \alpha[k] * Fk(x; \lambda[k])
#' }
#' 
#' where \eqn{\sum_i \alpha_i = 1}{sum(\alpha[i]) == 1}.
#'
#' @examples 
#' 
#' x <- rmixpois(1e5, c(5, 12, 19), c(1/3, 1/3, 1/3))
#' xx <- seq(-1, 50)
#' plot(prop.table(table(x)))
#' lines(xx, dmixpois(xx, c(5, 12, 19), c(1/3, 1/3, 1/3)), col = "red")
#' hist(pmixpois(x, c(5, 12, 19), c(1/3, 1/3, 1/3)))
#' 
#' xx <- seq(0, 50, by = 0.01)
#' plot(ecdf(x))
#' lines(xx, pmixpois(xx, c(5, 12, 19), c(1/3, 1/3, 1/3)), col = "red", lwd = 2)
#'
#' @name PoissonMix
#' @aliases PoissonMix
#' @aliases dmixpois
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Discrete
#'
#' @export

dmixpois <- function(x, lambda, alpha, log = FALSE) {
  
  if (is.vector(lambda))
    lambda <- matrix(lambda, nrow = 1)
  else if (!is.matrix(lambda))
    lambda <- as.matrix(lambda)
  
  if (is.vector(alpha))
    alpha <- matrix(alpha, nrow = 1)
  else if (!is.matrix(alpha))
    alpha <- as.matrix(alpha)
  
  cpp_dmixpois(x, lambda, alpha, log[1L])
}


#' @rdname PoissonMix
#' @export

pmixpois <- function(q, lambda, alpha, lower.tail = TRUE, log.p = FALSE) {
  
  if (is.vector(lambda))
    lambda <- matrix(lambda, nrow = 1)
  else if (!is.matrix(lambda))
    lambda <- as.matrix(lambda)
  
  if (is.vector(alpha))
    alpha <- matrix(alpha, nrow = 1)
  else if (!is.matrix(alpha))
    alpha <- as.matrix(alpha)
  
  cpp_pmixpois(q, lambda, alpha, lower.tail[1L], log.p[1L])
}


#' @rdname PoissonMix
#' @export

rmixpois <- function(n, lambda, alpha) {
  if (length(n) > 1) n <- length(n)
  
  if (is.vector(lambda))
    lambda <- matrix(lambda, nrow = 1)
  else if (!is.matrix(lambda))
    lambda <- as.matrix(lambda)
  
  if (is.vector(alpha))
    alpha <- matrix(alpha, nrow = 1)
  else if (!is.matrix(alpha))
    alpha <- as.matrix(alpha)
  
  cpp_rmixpois(n, lambda, alpha)
}

