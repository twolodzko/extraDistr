

#' Multinomial distribution
#'
#' Probability mass function and random generation
#' for the multinomial distribution.
#'
#' @param x 	 \eqn{k}-column matrix of quantiles.
#' @param n   number of observations. If \code{length(n) > 1},
#'             the length is taken to be the number required.
#' @param size numeric vector; number of trials (zero or more).
#' @param prob \eqn{k}-column numeric matrix; probability of success on each trial.
#' @param log  logical; if TRUE, probabilities p are given as log(p).
#'
#' @details
#'
#' Probability mass function
#' \deqn{
#' f(x) = \frac{n!}{\prod_{i=1}^k x_i} \prod_{i=1}^k p_i^{x_i}
#' }{
#' f(x) = n!/prod(x[i]!) * prod(p[i]^x[i])
#' }
#' 
#' @references 
#' Gentle, J.E. (2006). Random number generation and Monte Carlo methods. Springer.
#'
#' @seealso \code{\link[stats]{Binomial}}, \code{\link[stats]{Multinomial}}
#' 
#' @examples 
#' 
#' # Generating 10 random draws from multinomial distribution
#' # parametrized using a vector
#' 
#' (x <- rmnom(10, 3, c(1/3, 1/3, 1/3)))
#' 
#' # Results are consistent with dmultinom() from stats:
#' 
#' all.equal(dmultinom(x[1,], 3, c(1/3, 1/3, 1/3)),
#'           dmnom(x[1, , drop = FALSE], 3, c(1/3, 1/3, 1/3)))
#'
#' @name Multinomial
#' @aliases Multinomial
#' @aliases dmnom
#' 
#' @keywords distribution
#' @concept Multivariate
#' @concept Discrete
#'
#' @export

dmnom <- function(x, size, prob, log = FALSE) {
  if (is.vector(prob))
    prob <- matrix(prob, nrow = 1)
  else if (!is.matrix(prob))
    prob <- as.matrix(prob)
  
  if (is.vector(x))
    x <- matrix(x, nrow = 1)
  else if (!is.matrix(x))
    x <- as.matrix(x)
  
  cpp_dmnom(x, size, prob, log[1L])
}


#' @rdname Multinomial
#' @export

rmnom <- function(n, size, prob) {
  if (length(n) > 1)
    n <- length(n)
  
  if (is.vector(prob))
    prob <- matrix(prob, nrow = 1)
  else if (!is.matrix(prob))
    prob <- as.matrix(prob)
  
  cpp_rmnom(n, size, prob)
}

