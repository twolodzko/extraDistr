

#' Categorical distribution
#'
#' Probability mass function, distribution function, quantile function and random generation
#' for the categorical distribution.
#'
#' @param x,q	            vector of quantiles.
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
#' @examples 
#' 
#' # Generating 10 random draws from categorical distribution
#' # with k=3 categories occuring with equal probabilities
#' # parametrized using a vector
#' 
#' rcat(10, c(1/3, 1/3, 1/3))
#' 
#' # or with k=5 categories parametrized using a matrix of probabilities
#' # (generated from Dirichlet distribution)
#' 
#' p <- rdirichlet(10, c(1, 1, 1, 1, 1))
#' rcat(10, p)
#' 
#' x <- rcat(1e5, c(0.2, 0.4, 0.3, 0.1))
#' plot(prop.table(table(x)), type = "h")
#' lines(0:5, dcat(0:5, c(0.2, 0.4, 0.3, 0.1)), col = "red")
#'
#' @name Categorical
#' @aliases Categorical
#' @aliases dcat
#' @keywords distribution
#'
#' @export

dcat <- function(x, prob, log = FALSE) {
  if (is.vector(prob))
    prob <- matrix(prob, nrow = 1)
  else if (!is.matrix(prob))
    prob <- as.matrix(prob)
  .Call('extraDistr_cpp_dcat', PACKAGE = 'extraDistr', as.numeric(x), prob, log)
}


#' @rdname Categorical
#' @export

pcat <- function(q, prob, lower.tail = TRUE, log.p = FALSE) {
  if (is.vector(prob))
    prob <- matrix(prob, nrow = 1)
  else if (!is.matrix(prob))
    prob <- as.matrix(prob)
  .Call('extraDistr_cpp_pcat', PACKAGE = 'extraDistr', as.numeric(q), prob,
        lower.tail, log.p)
}


#' @rdname Categorical
#' @export

qcat <- function(p, prob, lower.tail = TRUE, log.p = FALSE, labels) {
  if (is.vector(prob))
    prob <- matrix(prob, nrow = 1)
  else if (!is.matrix(prob))
    prob <- as.matrix(prob)
  
  x <- .Call('extraDistr_cpp_qcat', PACKAGE = 'extraDistr', p, prob,
             lower.tail, log.p)
  
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
  
  if (is.vector(prob)) {
    k <- length(prob)
    if (any(prob < 0)) {
      warning("NaNs produced")
      x <- rep(NaN, n)
    } else {
      x <- sample.int(length(prob), size = n, replace = TRUE, prob = prob)
    }
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

