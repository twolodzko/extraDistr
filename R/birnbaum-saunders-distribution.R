

#' Birnbaum-Saunders (fatigue life) distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Birnbaum-Saunders (fatigue life) distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param alpha,beta,mu   shape, scale and location parameters.
#'                        Scale and shape must be positive.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'                        
#' @details
#' 
#' Probability density function
#' \deqn{
#' f(x) = \left (\frac{\sqrt{\frac{x-\mu} {\beta}} + \sqrt{\frac{\beta} 
#' {x-\mu}}} {2\alpha (x-\mu)} \right) 
#' \phi \left( \frac{1}{\alpha}\left( \sqrt{\frac{x-\mu}{\beta}} -
#' \sqrt{\frac{\beta}{x-\mu}} \right) \right)
#' }{
#' f(x) = ((sqrt((x-\mu)/\beta) + sqrt(\beta/(x-\mu)))/(2*\alpha*(x-\mu))) *
#' \phi((sqrt((x-\mu)/\beta) - sqrt(\beta/(x-\mu)))/\alpha)
#' }
#' 
#' Cumulative distribution function
#' \deqn{
#' F(x) = \Phi \left(\frac{1}{\alpha}\left( \sqrt{\frac{x-\mu}{\beta}} -
#' \sqrt{\frac{\beta}{x-\mu}} \right) \right)
#' }{
#' F(x) = \Phi(((sqrt((x-\mu)/\beta) - sqrt(\beta/(x-\mu)))/\alpha)
#' }
#' 
#' Quantile function
#' \deqn{
#' F^{-1}(p) = \left[\frac{\alpha}{2} \Phi^{-1}(p) +
#' \sqrt{\left(\frac{\alpha}{2} \Phi^{-1}(p)\right)^{2} + 1}\right]^{2} \beta + \mu
#' }{
#' F^-1(p) = (\alpha/2 * \Phi^-1(p) +
#' sqrt((\alpha/2 * \Phi^-1(p))^2 + 1)^2 * \beta + \mu
#' }
#'
#' @references
#' Birnbaum, Z. W. and Saunders, S. C. (1969).
#' A new family of life distributions.
#' Journal of Applied Probability, 6(2), 637-652.
#' 
#' @references 
#' Desmond, A. (1985) Stochastic models of failure in random environments.
#' Canadian Journal of Statistics, 13, 171-183.
#' 
#' @references 
#' Vilca-Labra, F., and Leiva-Sanchez, V. (2006).
#' A new fatigue life model based on the family of skew-elliptical distributions.
#' Communications in Statistics-Theory and Methods, 35(2), 229-244.
#' 
#' @references 
#' Leiva, V., Sanhueza, A., Sen, P. K., and Paula, G. A. (2008).
#' Random number generators for the generalized Birnbaum-Saunders distribution.
#' Journal of Statistical Computation and Simulation, 78(11), 1105-1118.
#'
#' @examples 
#' 
#' x <- rfatigue(1e5, .5, 2, 5)
#' hist(x, 100, freq = FALSE)
#' curve(dfatigue(x, .5, 2, 5), 2, 20, col = "red", add = TRUE)
#' hist(pfatigue(x, .5, 2, 5))
#' plot(ecdf(x))
#' curve(pfatigue(x, .5, 2, 5), 2, 20, col = "red", lwd = 2, add = TRUE)
#'
#' @name BirnbaumSaunders
#' @aliases BirnbaumSaunders
#' @aliases dfatigue
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dfatigue <- function(x, alpha, beta = 1, mu = 0, log = FALSE) {
  cpp_dfatigue(x, alpha, beta, mu, log[1L])
}


#' @rdname BirnbaumSaunders
#' @export

pfatigue <- function(q, alpha, beta = 1, mu = 0, lower.tail = TRUE, log.p = FALSE) {
  cpp_pfatigue(q, alpha, beta, mu, lower.tail[1L], log.p[1L])
}


#' @rdname BirnbaumSaunders
#' @export

qfatigue <- function(p, alpha, beta = 1, mu = 0, lower.tail = TRUE, log.p = FALSE) {
  cpp_qfatigue(p, alpha, beta, mu, lower.tail[1L], log.p[1L])
}


#' @rdname BirnbaumSaunders
#' @export

rfatigue <- function(n, alpha, beta = 1, mu = 0) {
  if (length(n) > 1) n <- length(n)
  cpp_rfatigue(n, alpha, beta, mu)
}

