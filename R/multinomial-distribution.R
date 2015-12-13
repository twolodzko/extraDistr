

#' Multinomial distribution
#'
#' Probability mass function and random generation
#' for the multinomial distribution.
#'
#' @param x 	 \eqn{m}-column matrix of quantiles.
#' @param n   number of observations. If \code{length(n) > 1},
#'             the length is taken to be the number required.
#' @param size numeric vector; number of trials (zero or more).
#' @param prob k-column numeric matrix; probability of success on each trial.
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
#' @seealso \code{\link{Binomial}}
#'
#' @name Multinomial
#' @aliases Multinomial
#' @aliases dmnom
#' @keywords distribution
#'
#' @export

dmnom <- function(x, size, prob, log = FALSE) {
  if (!is.matrix(prob))
    prob <- matrix(prob, nrow = 1)
  .Call('extraDistr_cpp_dmnom', PACKAGE = 'extraDistr', x, size, prob, log)
}


#' @rdname Multinomial
#' @export

rmnom <- function(n, size, prob) {
  if (length(n) > 1)
    n <- length(n)
  if (!is.matrix(prob))
    prob <- matrix(prob, nrow = 1)
  .Call('extraDistr_cpp_rmnom', PACKAGE = 'extraDistr', n, size, prob)
}

