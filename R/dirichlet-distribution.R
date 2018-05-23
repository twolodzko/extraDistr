

#' Dirichlet distribution
#'
#' Density function, cumulative distribution function and random generation
#' for the Dirichlet distribution.
#'
#' @param x               \eqn{k}-column matrix of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param alpha           \eqn{k}-values vector or \eqn{k}-column matrix;
#'                        concentration parameter. Must be positive.
#' @param log     	      logical; if TRUE, probabilities p are given as log(p).
#'
#' @details
#'
#' Probability density function
#' \deqn{
#' f(x) = \frac{\Gamma(\sum_k \alpha_k)}{\prod_k \Gamma(\alpha_k)} \prod_k x_k^{k-1}
#' }{
#' f(x) = \Gamma(sum(\alpha[k])) / prod(\Gamma(\alpha[k])) * prod(x[k]^{k-1})
#' }
#'
#' @references
#' Devroye, L. (1986). Non-Uniform Random Variate Generation. Springer-Verlag.
#' 
#' @examples
#' 
#' # Generating 10 random draws from Dirichlet distribution
#' # parametrized using a vector
#' 
#' rdirichlet(10, c(1, 1, 1, 1))
#' 
#' # or parametrized using a matrix where each row
#' # is a vector of parameters
#' 
#' alpha <- matrix(c(1, 1, 1, 1:3, 7:9), ncol = 3, byrow = TRUE)
#' rdirichlet(10, alpha)
#' 
#' @name Dirichlet
#' @aliases Dirichlet
#' @aliases ddirichlet
#' 
#' @keywords distribution
#' @concept Multivariate
#' @concept Continuous
#' 
#' @export

ddirichlet <- function(x, alpha, log = FALSE) {
  if (is.vector(alpha))
    alpha <- matrix(alpha, nrow = 1)
  else if (!is.matrix(alpha))
    alpha <- as.matrix(alpha)
  if (is.data.frame(x))
    x <- as.matrix(x)
  else if (is.vector(x))
    x <- matrix(x, byrow = TRUE, nrow = 1)
  cpp_ddirichlet(x, alpha, log[1L])
}


#' @rdname Dirichlet
#' @export

rdirichlet <- function (n, alpha) {
  if (length(n) > 1) n <- length(n)
  if (is.vector(alpha))
    alpha <- matrix(alpha, nrow = 1)
  else if (!is.matrix(alpha))
    alpha <- as.matrix(alpha)
  cpp_rdirichlet(n, alpha)
}

