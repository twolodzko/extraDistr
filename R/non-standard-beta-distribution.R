
#' Non-standard beta distribution
#' 
#' Non-standard form of beta distribution with lower and upper bounds
#' denoted as \code{min} and \code{max}. By default \code{min=0}
#' and \code{max=1} what leads to standard beta distribution.
#' 
#' @param x,q	          vector of quantiles.
#' @param p             vector of probabilities.
#' @param n             number of observations. If \code{length(n) > 1},
#'                      the length is taken to be the number required.
#' @param shape1,shape2 non-negative parameters of the Beta distribution.
#' @param min,max       lower and upper bounds.
#' @param log,log.p	    logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	  logical; if TRUE (default), probabilities are \eqn{P[X \leq x]},
#'                      otherwise, \eqn{P[X > x]}.
#'                      
#' @seealso \code{\link[stats]{Beta}}
#' 
#' @examples 
#' 
#' x <- rnsbeta(1e5, 5, 13, -4, 8)
#' hist(x, 100, freq = FALSE)
#' curve(dnsbeta(x, 5, 13, -4, 8), -4, 6, col = "red", add = TRUE) 
#' hist(pnsbeta(x, 5, 13, -4, 8))
#' plot(ecdf(x))
#' curve(pnsbeta(x, 5, 13, -4, 8), -4, 6, col = "red", lwd = 2, add = TRUE)
#' 
#' @name NSBeta
#' @aliases NSBeta
#' @aliases dgbeta
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'                     
#' @export

dnsbeta <- function(x, shape1, shape2, min = 0, max = 1, log = FALSE) {
  cpp_dnsbeta(x, shape1, shape2, min, max, log[1L])
}


#' @rdname NSBeta
#' @export

pnsbeta <- function(q, shape1, shape2, min = 0, max = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_pnsbeta(q, shape1, shape2, min, max, lower.tail[1L], log.p[1L])
}


#' @rdname NSBeta
#' @export

qnsbeta <- function(p, shape1, shape2, min = 0, max = 1, lower.tail = TRUE, log.p = FALSE) {
  cpp_qnsbeta(p, shape1, shape2, min, max, lower.tail[1L], log.p[1L])
}


#' @rdname NSBeta
#' @export

rnsbeta <- function(n, shape1, shape2, min = 0, max = 1) {
  if (length(n) > 1) n <- length(n)
  cpp_rnsbeta(n, shape1, shape2, min, max)
}

