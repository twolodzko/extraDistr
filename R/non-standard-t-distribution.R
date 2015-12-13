

#' Non-standard t-distribution
#' 
#' Probability mass function, distribution function and random generation
#' for non-standard t-distribution.
#'
#' @param x 	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param mu              vector of locations
#' @param sigma           vector of scale parameters.
#' @param df              degrees of freedom (> 0, maybe non-integer). \code{df = Inf} is allowed.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#' 
#' @name NonStandardT
#' @aliases NonStandardT
#' @aliases dnst
#' @keywords distribution
#' 
#' @export

dnst <- function(x, df, mu = 0, sigma = 1, log = FALSE) {
  .Call('extraDistr_cpp_dnst', PACKAGE = 'extraDistr', x, df, mu, sigma, log)
}


#' @rdname NonStandardT
#' @export

pnst <- function(x, df, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_pnst', PACKAGE = 'extraDistr', x, df, mu, sigma, lower.tail, log.p)
}


#' @rdname NonStandardT
#' @export

qnst <- function(p, df, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_qnst', PACKAGE = 'extraDistr', p, df, mu, sigma, lower.tail, log.p)
}


#' @rdname NonStandardT
#' @export

rnst <- function(n, df, mu = 0, sigma = 1) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rnst', PACKAGE = 'extraDistr', n, df, mu, sigma)
}




# dlst <- function(x, df, mu = 0, sigma = 1, log = FALSE) {
#   p <- dt((x - mu)/sigma, df)/sigma
#   if (log) log(p) else p
# }
# 
# plst <- function(q, df, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
#   pt((q - mu)/sigma, df, lower.tail = lower.tail, log.p = log.p)
# }
# 
# 
# qlst <- function(p, df, mu = 0, sigma = 1, lower.tail = TRUE, log.p = FALSE) {
#   qt(p, df, lower.tail = lower.tail, log.p = log.p)*sigma + mu
# }
# 
# rlst <- function(n, df, mu = 0, sigma = 1) {
#   rt(n, df)*sigma + mu
# }
