

#' Dirichlet-multinomial (multivariate \enc{Pólya}{Polya}) distribution
#'
#' Density function, cumulative distribution function and random generation
#' for the Dirichlet-multinomial (multivariate \enc{Pólya}{Polya}) distribution.
#'
#' @param x               \eqn{k}-column matrix of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param size            numeric vector; number of trials (zero or more).
#' @param alpha           \eqn{k}-values vector or \eqn{k}-column matrix;
#'                        concentration parameter. Must be positive.
#' @param log     	      logical; if TRUE, probabilities p are given as log(p).
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = \frac{\left(n!\right)\Gamma\left(\sum \alpha_k\right)}{\Gamma\left(n+\sum \alpha_k\right)}\prod_{k=1}^K\frac{\Gamma(x_{k}+\alpha_{k})}{\left(x_{k}!\right)\Gamma(\alpha_{k})}
#' }{
#' f(x) = (n! * \Gamma(sum(\alpha[k]))) / (\Gamma(n + sum(\alpha[k]))) * prod((\Gamma(x[k] + \alpha[k])) / (x[k]! * \Gamma(\alpha[k]))
#' }
#' 
#' @seealso \code{\link{Dirichlet}}, \code{\link{Multinomial}}
#'
#' @references 
#' Gentle, J.E. (2006). Random number generation and Monte Carlo methods. Springer.
#'
#' @references
#' Kvam, P. and Day, D. (2001) The multivariate Polya distribution in combat modeling.
#' Naval Research Logistics, 48, 1-17.
#' 
#' @name DirMnom
#' @aliases DirMnom
#' @aliases ddirmnom
#' @keywords distribution
#' 
#' @export

ddirmnom <- function(x, size, alpha, log = FALSE) {
  if (is.vector(alpha))
    alpha <- matrix(alpha, nrow = 1)
  else if (!is.matrix(alpha))
    alpha <- as.matrix(alpha)
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (is.vector(x))
    x <- matrix(x, byrow = TRUE, nrow = 1)
  .Call('extraDistr_cpp_ddirmnom', PACKAGE = 'extraDistr', x, size, alpha, log)
}


#' @rdname DirMnom
#' @export

rdirmnom <- function (n, size, alpha) {
  if (length(n) > 1) n <- length(n)
  if (is.vector(alpha))
    alpha <- matrix(alpha, nrow = 1)
  else if (!is.matrix(alpha))
    alpha <- as.matrix(alpha)
  .Call('extraDistr_cpp_rdirmnom', PACKAGE = 'extraDistr', n, size, alpha)
}

