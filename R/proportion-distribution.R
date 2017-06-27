

#' Beta distribution of proportions
#' 
#' Probability mass function, distribution function and random generation
#' for the reparametrized beta distribution.
#' 
#' @param x,q	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param size            non-negative real number; precision or number of binomial trials.
#' @param mean            mean proportion or probability of success on each trial;
#'                        \code{0 < mean < 1}.
#' @param prior           (see below) with \code{prior = 0} (default)
#'                        the distribution corresponds to re-parametrized beta
#'                        distribution used in beta regression. This parameter needs
#'                        to be non-negative.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#'                        
#' @details
#' 
#' Beta can be understood as a distribution of \eqn{x = k/\phi} proportions in
#' \eqn{\phi} trials where the average proportion is denoted as \eqn{\mu},
#' so it's parameters become \eqn{\alpha = \phi\mu} and
#' \eqn{\beta = \phi(1-\mu)} and it's density function becomes
#' 
#' \deqn{
#' f(x) = \frac{x^{\phi\mu+\pi-1} (1-x)^{\phi(1-\mu)+\pi-1}}{\mathrm{B}(\phi\mu+\pi, \phi(1-\mu)+\pi)}
#' }{
#' f(x) = (x^(\phi\mu+\pi-1) * (1-x)^(\phi(1-\mu)+\pi-1))/B(\phi\mu+\pi, \phi(1-\mu)+\pi)
#' }
#' 
#' where \eqn{\pi} is a \emph{prior} parameter, so the distribution is a
#' \emph{posterior} distribution after observing \eqn{\phi\mu} successes and
#' \eqn{\phi(1-\mu)} failures in \eqn{\phi} trials with binomial likelihood
#' and symmetric \eqn{\mathrm{Beta}(\pi, \pi)}{Beta(\pi, \pi)} prior for
#' probability of success. Parameter value \eqn{\pi = 1} corresponds to
#' uniform prior; \eqn{\pi = 1/2} corresponds to Jeffreys prior; \eqn{\pi = 0}
#' corresponds to "uninformative" Haldane prior, this is also the re-parametrized
#' distribution used in beta regression. With \eqn{\pi = 0} the distribution
#' can be understood as a continuous analog to binomial distribution dealing
#' with proportions rather then counts. Alternatively \eqn{\phi} may be
#' understood as precision parameter (as in beta regression).
#' 
#' Notice that in pre-1.8.4 versions of this package, \code{prior} was not settable
#' and by default fixed to one, instead of zero. To obtain the same results as in
#' the previous versions, use \code{prior = 1} in each of the functions. 
#' 
#' @seealso
#' \code{\link{beta}}, \code{\link{binomial}}
#' 
#' @references
#' Ferrari, S., & Cribari-Neto, F. (2004). Beta regression for modelling rates and proportions.
#' Journal of Applied Statistics, 31(7), 799-815.
#' 
#' @references 
#' Smithson, M., & Verkuilen, J. (2006). A better lemon squeezer? Maximum-likelihood regression
#' with beta-distributed dependent variables.
#' Psychological Methods, 11(1), 54-71.
#' 
#' @examples 
#' 
#' x <- rprop(1e5, 100, 0.33)
#' hist(x, 100, freq = FALSE)
#' curve(dprop(x, 100, 0.33), 0, 1, col = "red", add = TRUE)
#' hist(pprop(x, 100, 0.33))
#' plot(ecdf(x))
#' curve(pprop(x, 100, 0.33), 0, 1, col = "red", lwd = 2, add = TRUE)
#' 
#' n <- 500
#' p <- 0.23
#' k <- rbinom(1e5, n, p)
#' hist(k/n, freq = FALSE, 100)
#' curve(dprop(x, n, p), 0, 1, col = "red", add = TRUE, n = 500)
#'                        
#' @name PropBeta
#' @aliases PropBeta
#' @aliases dprop
#' 
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#'
#' @export

dprop <- function(x, size, mean, prior = 0, log = FALSE) {
  cpp_dprop(x, size, mean, prior, log[1L])
}


#' @rdname PropBeta
#' @export

pprop <- function(q, size, mean, prior = 0, lower.tail = TRUE, log.p = FALSE) {
  cpp_pprop(q, size, mean, prior, lower.tail[1L], log.p[1L])
}


#' @rdname PropBeta
#' @export

qprop <- function(p, size, mean, prior = 0, lower.tail = TRUE, log.p = FALSE) {
  cpp_qprop(p, size, mean, prior, lower.tail[1L], log.p[1L])
}


#' @rdname PropBeta
#' @export

rprop <- function(n, size, mean, prior = 0) {
  if (length(n) > 1) n <- length(n)
  cpp_rprop(n, size, mean, prior)
}

