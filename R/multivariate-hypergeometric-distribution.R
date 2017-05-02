

#' Multivariate hypergeometric distribution
#'
#' Probability mass function and random generation
#' for the multivariate hypergeometric distribution.
#'
#' @param x 	 \eqn{m}-column matrix of quantiles.
#' @param nn   number of observations. If \code{length(n) > 1},
#'             the length is taken to be the number required.
#' @param n    \eqn{m}-length vector or \eqn{m}-column matrix
#'             of numbers of balls in \eqn{m} colors.
#' @param k    the number of balls drawn from the urn.
#' @param log  logical; if TRUE, probabilities p are given as log(p).
#'
#' @details
#'
#' Probability mass function
#' \deqn{
#' f(x) = \frac{\prod_{i=1}^m {n_i \choose x_i}}{{N \choose k}}
#' }{
#' f(x) = prod(choose(n, x)) / choose(N, k)
#' }
#' 
#' @details
#' The multivariate hypergeometric distribution is generalization of
#' hypergeometric distribution. It is used for sampling \emph{without} replacement
#' \eqn{k} out of \eqn{N} marbles in \eqn{m} colors, where each of the colors appears
#' \eqn{n_i}{n[i]} times. Where \eqn{k=\sum_{i=1}^m x_i}{k=sum(x)},
#' \eqn{N=\sum_{i=1}^m n_i}{N=sum(n)} and \eqn{k \le N}{k<=N}.
#' 
#' @references 
#' Gentle, J.E. (2006). Random number generation and Monte Carlo methods. Springer.
#'
#' @seealso \code{\link[stats]{Hypergeometric}}
#' 
#' @examples 
#' 
#' # Generating 10 random draws from multivariate hypergeometric
#' # distribution parametrized using a vector
#' 
#' rmvhyper(10, c(10, 12, 5, 8, 11), 33)
#'
#' @name MultiHypergeometric
#' @aliases MultiHypergeometric
#' @aliases dmvhyper
#' 
#' @keywords distribution
#' @concept Multivariate
#' @concept Discrete
#'
#' @export

dmvhyper <- function(x, n, k, log = FALSE) {
  if (is.vector(n))
    n <- matrix(n, nrow = 1)
  else if (!is.matrix(n))
    n <- as.matrix(n)
  
  if (is.vector(x))
    x <- matrix(x, nrow = 1)
  else if (!is.matrix(x))
    x <- as.matrix(x)
  
  cpp_dmvhyper(x, n, k, log[1L])
}


#' @rdname MultiHypergeometric
#' @export

rmvhyper <- function(nn, n, k) {
  if (length(nn) > 1)
    nn <- length(nn)
  
  if (is.vector(n) && length(k) == 1) {
    if (anyNA(n) || is.na(k)) {
      warning("NAs produced")
      return(matrix(rep(NA, nn), nrow = nn, byrow = TRUE))
    }
    if (sum(n) == k)
      return(matrix(rep(n, nn), nrow = nn, byrow = TRUE))
    n <- matrix(n, nrow = 1)
  } else if (!is.matrix(n)) {
    n <- as.matrix(n)
  }
  
  cpp_rmvhyper(nn, n, k)
}

