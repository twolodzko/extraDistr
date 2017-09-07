

#' Bivariate Poisson distribution
#'
#' Probability mass function and random generation for the bivariate Poisson distribution.
#'
#' @param x,y	  vectors of quantiles; alternatively x may be a two-column
#'              matrix (or data.frame) and y may be omitted.
#' @param n	    number of observations. If \code{length(n) > 1},
#'              the length is taken to be the number required.
#' @param a,b,c positive valued parameters.
#' @param log   logical; if TRUE, probabilities p are given as log(p).
#'
#' @details
#'
#' Probability mass function
#' \deqn{
#' f(x) = \exp \{-(a+b+c)\} \frac{a^x}{x!} \frac{b^y}{y!} \sum_{k=0}^{\min(x,y)}
#' {x \choose k} {y \choose k} k! \left( \frac{c}{ab} \right)^k
#' }{
#' f(x) = exp(-(a+b+c)) * (a^x)/x! * (b^y)/y! *
#' sum(choose(x,k)*choose(y,k)*k!*(c/(a*b))^k)
#' }
#' 
#' @references 
#' Karlis, D. and Ntzoufras, I. (2003). Analysis of sports data by using bivariate Poisson models.
#' Journal of the Royal Statistical Society: Series D (The Statistician), 52(3), 381-393.
#' 
#' @references 
#' Kocherlakota, S. and Kocherlakota, K. (1992) Bivariate Discrete Distributions.
#' New York: Dekker.
#' 
#' @references 
#' Johnson, N., Kotz, S. and Balakrishnan, N. (1997). Discrete Multivariate Distributions.
#' New York: Wiley.
#' 
#' @references 
#' Holgate, P. (1964). Estimation for the bivariate Poisson distribution.
#' Biometrika, 51(1-2), 241-287.
#' 
#' @references 
#' Kawamura, K. (1984). Direct calculation of maximum likelihood estimator for the bivariate
#' Poisson distribution. Kodai mathematical journal, 7(2), 211-221.
#' 
#' @examples 
#' 
#' x <- rbvpois(5000, 7, 8, 5)
#' image(prop.table(table(x[,1], x[,2])))
#' colMeans(x)
#'
#' @seealso \code{\link[stats]{Poisson}}
#'
#' @name BivPoiss
#' @aliases BivPoiss
#' @aliases dbvpois
#' 
#' @keywords distribution
#' @concept Multivariate
#' @concept Discrete
#'
#' @export

dbvpois <- function(x, y = NULL, a, b, c, log = FALSE) {
  if (is.null(y)) {
    if ((is.matrix(x) || is.data.frame(x)) && ncol(x) == 2) {
      y <- x[, 2]
      x <- x[, 1]
    } else {
      stop("y is not provided while x is not a two-column matrix")
    }
  }
  cpp_dbpois(x, y, a, b, c, log[1L])
}


#' @rdname BivPoiss
#' @export

rbvpois <- function(n, a, b, c) {
  if (length(n) > 1) n <- length(n)
  cpp_rbpois(n, a, b, c)
}

