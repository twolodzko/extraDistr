

#' Negative hypergeometric distribution
#'
#' Probability mass function, distribution function, quantile function and random generation
#' for the negative hypergeometric distribution.
#'
#' @param x,q        vector of quantiles representing the number of balls drawn without
#'                   replacement from an urn which contains both black and white balls.
#' @param p	         vector of probabilities.
#' @param nn         number of observations. If \code{length(n) > 1},
#'                   the length is taken to be the number required.
#' @param m          the number of white balls in the urn.
#' @param n          the number of black balls in the urn.
#' @param r          the number of white balls that needs to be drawn for the sampling
#'                   to be stopped.
#' @param log,log.p  logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                   otherwise, \eqn{P[X > x]}.
#'                        
#'
#' @details
#'
#' Negative hypergeometric distribution describes number of balls \eqn{x} observed
#' until drawing without replacement to obtain \eqn{r} white balls
#' from the urn containing \eqn{m} white balls and \eqn{n} black balls,
#' and is defined as
#' 
#' \deqn{
#' f(x) = \frac{{x-1 \choose r-1}{m+n-x \choose m-r}}{{m+n \choose n}}
#' }{
#' f(x) = choose(x-1, r-1)*choose(m+n-x, m-r)/choose(m+n, n)
#' }
#' 
#' The algorithm used for calculating probability mass function,
#' cumulative distribution function and quantile function is based
#' on Fortran program NHYPERG created by Berry and Mielke (1996, 1998).
#' Random generation is done by inverse transform sampling.
#' 
#' @references 
#' Berry, K. J., & Mielke, P. W. (1998).
#' The negative hypergeometric probability distribution:
#' Sampling without replacement from a finite population.
#' Perceptual and motor skills, 86(1), 207-210.
#' \url{http://pms.sagepub.com/content/86/1/207.full.pdf}
#' 
#' @references 
#' Berry, K. J., & Mielke, P. W. (1996).
#' Exact confidence limits for population proportions based on the negative
#' hypergeometric probability distribution.
#' Perceptual and motor skills, 83(3 suppl), 1216-1218.
#' \url{http://pms.sagepub.com/content/83/3_suppl/1216.full.pdf}
#' 
#' @references 
#' Schuster, E. F., & Sype, W. R. (1987).
#' On the negative hypergeometric distribution.
#' International Journal of Mathematical Education in Science and Technology, 18(3), 453-459.
#' 
#' @references 
#' Chae, K. C. (1993).
#' Presenting the negative hypergeometric distribution to the introductory statistics courses.
#' International Journal of Mathematical Education in Science and Technology, 24(4), 523-526.
#' 
#' @references 
#' Jones, S.N. (2013). A Gaming Application of the Negative Hypergeometric Distribution.
#' UNLV Theses, Dissertations, Professional Papers, and Capstones. Paper 1846.
#' \url{http://digitalscholarship.unlv.edu/cgi/viewcontent.cgi?article=2847&context=thesesdissertations}
#'
#' @seealso \code{\link[stats]{Hypergeometric}}
#' 
#' @examples
#' 
#' x <- rnhyper(1e5, 60, 35, 15)
#' xx <- 15:95
#' plot(prop.table(table(x)))
#' lines(xx, dnhyper(xx, 60, 35, 15), col = "red")
#' hist(pnhyper(x, 60, 35, 15))
#' 
#' xx <- seq(0, 100, by = 0.01)
#' plot(ecdf(x))
#' lines(xx, pnhyper(xx, 60, 35, 15), col = "red", lwd = 2)
#'
#' @name NegHyper
#' @aliases NegHyper
#' @aliases dnhyper
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Discrete
#'
#' @export

dnhyper <- function(x, n, m, r, log = FALSE) {
  cpp_dnhyper(x, n, m, r, log[1L])
}


#' @rdname NegHyper
#' @export

pnhyper <- function(q, n, m, r, lower.tail = TRUE, log.p = FALSE) {
  cpp_pnhyper(q, n, m, r, lower.tail[1L], log.p[1L])
}


#' @rdname NegHyper
#' @export

qnhyper <- function(p, n, m, r, lower.tail = TRUE, log.p = FALSE) {
  cpp_qnhyper(p, n, m, r, lower.tail[1L], log.p[1L])
}


#' @rdname NegHyper
#' @export

rnhyper <- function(nn, n, m, r) {
  if (length(nn) > 1)
    nn <- length(nn)
  cpp_rnhyper(nn, n, m, r)
}

