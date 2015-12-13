
#' General beta distribution
#' 
#' General form of beta distribution with lower and upper bounds
#' denoted as \code{min} and \code{max}. By default \code{min=0}
#' and \code{max=1} what leads to standard beta distribution.
#' 
#' @param x             vector of quantiles.
#' @param p             vector of probabilities.
#' @param n             number of observations. If \code{length(n) > 1},
#'                      the length is taken to be the number required.
#' @param shape1,shape2 non-negative parameters of the Beta distribution.
#' @param min,max       lower and upper bounds.
#' @param log,log.p	    logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	  logical; if TRUE (default), probabilities are \eqn{P[X \leq x]},
#'                      otherwise, \eqn{P[X > x]}.
#'                      
#' @seealso \code{\link{Beta}}
#' 
#' @name GeneralBeta
#' @aliases GeneralBeta
#' @aliases dgbeta
#' @keywords distribution
#'                     
#' @export

dgbeta <- function(x, shape1, shape2, min = 0, max = 1, log = FALSE) {
  x <- (x-min)/(max-min)
  p <- dbeta(x, shape1, shape2, log = log)
  if (log) p-log(max-min)
  else     p/(max-min)
}


#' @rdname GeneralBeta
#' @export

pgbeta <- function(x, shape1, shape2, min = 0, max = 1,
                   lower.tail = TRUE, log.p = FALSE) {
  x <- (x-min)/(max-min)
  pbeta(x, shape1, shape2, lower.tail = lower.tail, log.p = log.p)
}


#' @rdname GeneralBeta
#' @export

qgbeta <- function(p, shape1, shape2, min = 0, max = 1,
                   lower.tail = TRUE, log.p = FALSE) {
  x <- qbeta(p, shape1, shape2, lower.tail = lower.tail, log.p = log.p)
  x * (max-min) + min
}


#' @rdname GeneralBeta
#' @export

rgbeta <- function(n, shape1, shape2, min = 0, max = 1, log = FALSE) {
  x <- rbeta(n, shape1, shape2)
  x * (max-min) + min
}

