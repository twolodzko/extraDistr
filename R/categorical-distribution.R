

#' Categorical distribution
#'
#' Probability mass function, distribution function, quantile function and random generation
#' for the categorical distribution.
#'
#' @param x 	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param prob            vector of length \eqn{k}, or \eqn{k}-column matrix
#'                        of probabilities. Probabilities need to sum up to 1.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#' @param labels          if provided, labeled \code{factor} vector is returned.
#'                        Number of labels needs to be the same as
#'                        number of categories (number of columns in prob).
#'
#' @name Categorical
#' @aliases Categorical
#' @aliases dcat
#' @keywords distribution
#'
#' @export

dcat <- function(x, prob, log = FALSE) {
  if (!is.matrix(prob))
    prob <- matrix(prob, nrow = 1)
  .Call('extraDistr_cpp_dcat', PACKAGE = 'extraDistr', x, prob, log)
}


#' @rdname Categorical
#' @export

pcat <- function(x, prob, lower.tail = TRUE, log.p = FALSE) {
  if (is.vector(prob))
    prob <- matrix(prob, nrow = 1)
  else if (!is.matrix(prob))
    prob <- as.matrix(prob)
  .Call('extraDistr_cpp_pcat', PACKAGE = 'extraDistr', x, prob, lower.tail, log.p)
}


#' @rdname Categorical
#' @export

qcat <- function(p, prob, lower.tail = TRUE, log.p = FALSE, labels) {
  if (is.vector(prob))
    prob <- matrix(prob, nrow = 1)
  else if (!is.matrix(prob))
    prob <- as.matrix(prob)
  
  x <- .Call('extraDistr_cpp_qcat', PACKAGE = 'extraDistr', p, prob, lower.tail, log.p)
  
  if (!missing(labels)) {
    if (length(labels) != ncol(prob))
      warning("Wrong number of labels.")
    else
      return(factor(x, levels = 1:ncol(prob), labels = labels))
  }
  
  return(x)
}


#' @rdname Categorical
#' @export

rcat <- function(n, prob, labels) {
  if (length(n) > 1) n <- length(n)
  
  if (!is.matrix(prob)) {
    k <- length(prob)
    x <- sample.int(length(prob), size = n, replace = TRUE, prob = prob)
  } else {
    k <- ncol(prob)
    x <- .Call('extraDistr_cpp_rcat', PACKAGE = 'extraDistr', n, prob)
  }
  
  if (!missing(labels)) {
    if (length(labels) != k)
      warning("Wrong number of labels.")
    else
      return(factor(x, levels = 1:k, labels = labels))
  }
  
  return(x)
}

