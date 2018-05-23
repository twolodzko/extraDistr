

#' Shifted Gompertz distribution
#'
#' Density, distribution function, and random generation
#' for the shifted Gompertz distribution.
#'
#' @param x,q	            vector of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param b,eta           positive valued scale and shape parameters;
#'                        both need to be positive.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'
#' @details
#' 
#' If \eqn{X} follows exponential distribution parametrized by scale \eqn{b} and
#' \eqn{Y} follows reparametrized Gumbel distribution with cumulative distribution function
#' \eqn{F(x) = \exp(-\eta  e^{-bx})}{F(x) = exp(-\eta*exp(-b*x))} parametrized by 
#' scale \eqn{b} and shape \eqn{\eta}, then \eqn{\max(X,Y)}{max(X,Y)} follows shifted
#' Gompertz distribution parametrized by scale \eqn{b>0} and shape \eqn{\eta>0}.
#' The above relation is used by \code{rsgomp} function for random generation from
#' shifted Gompertz distribution.
#'
#' Probability density function
#' \deqn{
#' f(x) = b e^{-bx} \exp(-\eta e^{-bx}) \left[1 + \eta(1 - e^{-bx})\right]
#' }{
#' f(x) = b*exp(-b*x) * exp(-\eta*exp(-b*x)) * (1 + \eta*(1 - exp(-b*x)))
#' }
#'
#' Cumulative distribution function
#' \deqn{
#' F(x) = (1-e^{-bx}) \exp(-\eta e^{-bx})
#' }{
#' F(x) = (1-exp(-b*x)) * exp(-\eta*exp(-b*x))
#' }
#' 
#' @references 
#' Bemmaor, A.C. (1994).
#' Modeling the Diffusion of New Durable Goods: Word-of-Mouth Effect Versus Consumer Heterogeneity.
#' [In:] G. Laurent, G.L. Lilien & B. Pras. Research Traditions in Marketing.
#' Boston: Kluwer Academic Publishers. pp. 201-223.
#' 
#' @references 
#' Jimenez, T.F. and Jodra, P. (2009).
#' A Note on the Moments and Computer Generation of the Shifted Gompertz Distribution.
#' Communications in Statistics - Theory and Methods, 38(1), 78-89.
#' 
#' @references
#' Jimenez T.F. (2014).
#' Estimation of the Parameters of the Shifted Gompertz Distribution,
#' Using Least Squares, Maximum Likelihood and Moments Methods.
#' Journal of Computational and Applied Mathematics, 255(1), 867-877.
#' 
#' @examples 
#' 
#' x <- rsgomp(1e5, 0.4, 1)
#' hist(x, 50, freq = FALSE)
#' curve(dsgomp(x, 0.4, 1), 0, 30, col = "red", add = TRUE)
#' hist(psgomp(x, 0.4, 1))
#' plot(ecdf(x))
#' curve(psgomp(x, 0.4, 1), 0, 30, col = "red", lwd = 2, add = TRUE)
#'
#' @name ShiftGomp
#' @aliases ShiftGomp
#' @aliases dsgomp
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dsgomp <- function(x, b, eta, log = FALSE) {
  cpp_dsgomp(x, b, eta, log[1L])
}


#' @rdname ShiftGomp
#' @export

psgomp <- function(q, b, eta, lower.tail = TRUE, log.p = FALSE) {
  cpp_psgomp(q, b, eta, lower.tail[1L], log.p[1L])
}


#' @rdname ShiftGomp
#' @export

rsgomp <- function(n, b, eta) {
  if (length(n) > 1) n <- length(n)
  cpp_rsgomp(n, b, eta)
}

