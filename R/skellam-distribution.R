

#' Skellam distribution
#'
#' Probability mass function and random generation
#' for the Skellam distribution.
#'
#' @param x 	            vector of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param mu1,mu2         positive valued parameters.
#' @param log     	      logical; if TRUE, probabilities p are given as log(p).
#'
#' @details
#' 
#' If \eqn{X} and \eqn{Y} follow Poisson distributions with means
#' \eqn{\mu_1}{\mu[1]} and \eqn{\mu_2}{\mu[2]}, than \eqn{X-Y} follows
#' Skellam distribution.
#'
#' Probability mass function
#' \deqn{
#' f(x) = e^{-(\mu_1\!+\!\mu_2)} \left(\frac{\mu_1}{\mu_2}\right)^{k/2}\!\!I_{k}(2\sqrt{\mu_1\mu_2})
#' }{
#' f(x) = exp(-(\mu1+\mu2)) * (\mu1/\mu2)^(x/2) * besselI(2*sqrt(\mu1*\mu2), x)
#' }
#' 
#' @references
#' Karlis, D., & Ntzoufras, I. (2006). Bayesian analysis of the differences of count data.
#' Statistics in medicine, 25(11), 1885-1905.
#' 
#' @examples 
#' 
#' x <- rskellam(1e5, 5, 13)
#' xx <- -40:40
#' plot(prop.table(table(x)), type = "h")
#' lines(xx, dskellam(xx, 5, 13), col = "red")
#'
#' @name Skellam
#' @aliases Skellam
#' @aliases dskellam
#' @keywords distribution
#'
#' @export

dskellam <- function(x, mu1, mu2, log = FALSE) {
  .Call('extraDistr_cpp_dskellam', PACKAGE = 'extraDistr', x, mu1, mu2, log)
}


#' @rdname Skellam
#' @export

rskellam <- function(n, mu1, mu2) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rskellam', PACKAGE = 'extraDistr', n, mu1, mu2)
}

