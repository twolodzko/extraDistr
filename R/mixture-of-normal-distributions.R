

#' Mixture of normal distributions
#'
#' Density, distribution function and random generation
#' for the mixture of normal distributions.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param mean 	          matrix (or vector) of means.
#' @param sd              matrix (or vector) of standard deviations.
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
#' f(x) = \alpha_1 f_1(x; \mu_1, \sigma_1) + \dots + \alpha_k f_k(x; \mu_k, \sigma_k)
#' }{
#' f(x) = \alpha[1] * f1(x; \mu[1], \sigma[1]) + \dots + \alpha[k] * fk(x; \mu[k], \sigma[k])
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = \alpha_1 F_1(x; \mu_1, \sigma_1) + \dots + \alpha_k F_k(x; \mu_k, \sigma_k)
#' }{
#' F(x) = \alpha[1] * F1(x; \mu[1], \sigma[1]) + \dots + \alpha[k] * Fk(x; \mu[k], \sigma[k])
#' }
#' 
#' where \eqn{\sum_i \alpha_i = 1}{sum(\alpha[i]) == 1}.
#'
#' @examples 
#' 
#' x <- rmixnorm(1e5, c(0.5, 3, 6), c(3, 1, 1), c(1/3, 1/3, 1/3))
#' hist(x, 100, freq = FALSE)
#' curve(dmixnorm(x, c(0.5, 3, 6), c(3, 1, 1), c(1/3, 1/3, 1/3)),
#'       -20, 20, n = 500, col = "red", add = TRUE)
#' hist(pmixnorm(x, c(0.5, 3, 6), c(3, 1, 1), c(1/3, 1/3, 1/3)))
#' plot(ecdf(x))
#' curve(pmixnorm(x, c(0.5, 3, 6), c(3, 1, 1), c(1/3, 1/3, 1/3)),
#'       -20, 20, n = 500, col = "red", lwd = 2, add = TRUE)
#'
#' @name NormalMix
#' @aliases NormalMix
#' @aliases dmixnorm
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dmixnorm <- function(x, mean, sd, alpha, log = FALSE) {
  
  if (is.vector(mean))
    mean <- matrix(mean, nrow = 1)
  else if (!is.matrix(mean))
    mean <- as.matrix(mean)
  
  if (is.vector(sd))
    sd <- matrix(sd, nrow = 1)
  else if (!is.matrix(sd))
    sd <- as.matrix(sd)
  
  if (is.vector(alpha))
    alpha <- matrix(alpha, nrow = 1)
  else if (!is.matrix(alpha))
    alpha <- as.matrix(alpha)
  
  cpp_dmixnorm(x, mean, sd, alpha, log[1L])
}


#' @rdname NormalMix
#' @export

pmixnorm <- function(q, mean, sd, alpha, lower.tail = TRUE, log.p = FALSE) {
  
  if (is.vector(mean))
    mean <- matrix(mean, nrow = 1)
  else if (!is.matrix(mean))
    mean <- as.matrix(mean)
  
  if (is.vector(sd))
    sd <- matrix(sd, nrow = 1)
  else if (!is.matrix(sd))
    sd <- as.matrix(sd)
  
  if (is.vector(alpha))
    alpha <- matrix(alpha, nrow = 1)
  else if (!is.matrix(alpha))
    alpha <- as.matrix(alpha)
  
  cpp_pmixnorm(q, mean, sd, alpha, lower.tail[1L], log.p[1L])
}


#' @rdname NormalMix
#' @export

rmixnorm <- function(n, mean, sd, alpha) {
  if (length(n) > 1) n <- length(n)
  
  if (is.vector(mean))
    mean <- matrix(mean, nrow = 1)
  else if (!is.matrix(mean))
    mean <- as.matrix(mean)
  
  if (is.vector(sd))
    sd <- matrix(sd, nrow = 1)
  else if (!is.matrix(sd))
    sd <- as.matrix(sd)
  
  if (is.vector(alpha))
    alpha <- matrix(alpha, nrow = 1)
  else if (!is.matrix(alpha))
    alpha <- as.matrix(alpha)
  
  cpp_rmixnorm(n, mean, sd, alpha)
}

