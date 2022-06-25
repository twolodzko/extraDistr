

#' Truncated negative binomial distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the truncated negative binomial distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param size	          target for number of successful trials, or dispersion
#'                        parameter (the shape parameter of the gamma mixing
#'                        distribution). Must be strictly positive, need not be
#'                        integer.
#' @param prob            probability of success in each trial. \code{0 < prob <= 1}.
#' @param mu              alternative parameterization via mean
#' @param a,b             lower and upper truncation points (\code{a < x <= b}).
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @references
#' Hilbe, J. (2011). Censored and truncated count models. In 
#' *Negative Binomial Regression* (pp. 387-406). Cambridge: Cambridge University 
#' Press. \url{https://doi.org/10.1017/CBO9780511973420.013}
#' 
#' @seealso \code{\link[stats]{NegBinomial}}
#' 
#' @examples 
#' 
#' # Right-truncated negative binomial
#' ## random sample
#' x <- rtnbinom(1e5, size = 2, prob = 0.1, b = 25)
#' plot(prop.table(table(x)))
#' 
#' ## distribution
#' xx <- seq(-1, 30)
#' lines(xx, dtnbinom(xx, size = 2, prob = 0.1, b = 25), col = "red")
#' 
#' hist(ptnbinom(x, size = 2, prob = 0.1, b = 25), breaks = 35)
#' 
#' xx <- seq(0, 30, by = 0.01)
#' plot(ecdf(x))
#' lines(xx, ptnbinom(xx, size = 2, prob = 0.1, b = 25), col = "red", lwd = 2)
#' 
#' uu <- seq(0, 1, by = 0.001)
#' lines(qtnbinom(uu, size = 2, prob = 0.1, b = 25), uu, col = "blue", lty = 2)
#' 
#' # Zero-truncated negative binomial (mu parameterization)
#' ## random sample
#' x <- rtnbinom(1e5, size = 2, mu = 5, a = 0)
#' plot(prop.table(table(x)))
#' 
#' ## distribution
#' xx <- seq(-1, 50)
#' lines(xx, dtnbinom(xx, size = 2, mu = 5, a = 0), col = "red")
#' hist(ptnbinom(x, size = 2, mu = 5, a = 0))
#' 
#' xx <- seq(0, 50, by = 0.01)
#' plot(ecdf(x))
#' lines(xx, ptnbinom(xx, size = 2, mu = 5, a = 0), col = "red", lwd = 2)
#' lines(qtnbinom(uu, size = 2, mu = 5, a = 0), uu, col = "blue", lty = 2)
#'
#' @name TruncNegBinom
#' @aliases TruncNegBinom
#' @aliases TruncNB
#' @aliases dtnbinom
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Discrete
#'
#' @export

dtnbinom <- function(x, size, prob, mu, a = -Inf, b = Inf, log = FALSE) {
  if (!missing(mu)) {
    if(!missing(prob))
      stop("'prob' and 'mu' both specified")
    cpp_dtnbinom_mu(x, size, mu, a, b, log)
  }
  else cpp_dtnbinom(x, size, prob, a, b, log)
}


#' @rdname TruncNegBinom
#' @export

ptnbinom <- function(q, size, prob, mu, a = -Inf, b = Inf, lower.tail = TRUE, log.p = FALSE) {
  if (!missing(mu)) {
    if(!missing(prob))
      stop("'prob' and 'mu' both specified")
    cpp_ptnbinom_mu(q, size, mu, a, b, lower.tail[1L], log.p[1L])
  }
  else cpp_ptnbinom(q, size, prob, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname TruncNegBinom
#' @export

qtnbinom <- function(p, size, prob, mu, a = -Inf, b = Inf, lower.tail = TRUE, log.p = FALSE) {
  if (!missing(mu)) {
    if (!missing(prob)) 
      stop("'prob' and 'mu' both specified")
    cpp_qtnbinom_mu(p, size, mu, a, b, lower.tail[1L], log.p[1L])
  }
  else cpp_qtnbinom(p, size, prob, a, b, lower.tail[1L], log.p[1L])
}


#' @rdname TruncNegBinom
#' @export

rtnbinom <- function(n, size, prob, mu, a = -Inf, b = Inf) {
  if (length(n) > 1) n <- length(n)
  if (!missing(mu)) {
    if (!missing(prob)) 
      stop("'prob' and 'mu' both specified")
    cpp_rtnbinom_mu(n, size, mu, a, b)
  }
  else cpp_rtnbinom(n, size, prob, a, b)
}
