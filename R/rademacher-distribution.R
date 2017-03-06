

#' Random generation from Rademacher distribution
#'
#' Random generation for the Rademacher distribution (values -1 and +1
#' with equal probability).
#'
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#'
#' @name Rademacher
#' @aliases Rademacher
#' @aliases rsign
#'
#' @keywords distribution
#' 
#' @export

rsign <- function(n) {
  if (length(n) > 1) n <- length(n)
  cpp_rsign(n)
}

