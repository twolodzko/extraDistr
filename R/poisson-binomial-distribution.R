

#' Poisson-binomial distribution
#'
#' Probability mass function, distribution function and random generation
#' for the Poisson-binomial distribution.
#'
#' @param x 	            vector of quantiles.
#' @param p	              vector of probabilities.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param prob            parameters.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#' @param nsim            number of samples in Monte Carlo simulation for calculating
#'                        density and distribution function; the higher is more precise
#'                        but slower.
#' @param method          method for approximating density and distribution function;
#'                        "MC" for Monte Carlo simulation, "PA" for poisson approximation,
#'                        "BA" for Binomial approximation, "NA" for Normal approximation,
#'                        "RNA" for refined Normal approximation.
#'
#' @details
#' There are multiple ways of estimating probability mass function and distribution function
#' of Poisson-Binomial, however direct methods are numerically unstable. Multiple indirect
#' and approximate methods were proposed (Chen and Liu, 1997, Hong, 2013, Volkowa, 1996,
#' Choi and Xia, 2002). This function enables to use Monte Carlo simulation and
#' approximation-based methods. Poisson-Binomial can be approximated using Poisson
#' distribution by Le Cam theorem, however Binomial approximation is generally better
#' (Chen and Liu, 2002). Normal approximation or its correction proposed by Volkova (1996)
#' are also possible. Other methods are available in \pkg{poibin} package by Yili Hong.
#'
#' @references
#' Chen, L.H.Y. (1974). On the Convergence of Poisson Binomial to Poisson
#' Distributions. The Annals of Probability, 2(1), 178-180.
#'
#' @references
#' Chen, S.X. and Liu, J.S. (1997). Statistical applications of the Poisson-Binomial
#' and conditional Bernoulli distributions. Statistica Sinica 7, 875-892.
#'
#' @references
#' Chen, S.X. (1993). Poisson-Binomial distribution, conditional Bernoulli distribution
#' and maximum entropy. Technical Report. Department of Statistics, Harvard University.
#'
#' @references
#' Chen, X.H., Dempster, A.P. and Liu, J.S. (1994). Weighted finite population
#' sampling to maximize entropy. Biometrika 81, 457-469.
#'
#' @references
#' Wang, Y.H. (1993). On the number of successes in independent trials.
#' Statistica Sinica 3(2): 295-312.
#'
#' @references
#' Hong, Y. (2013). On computing the distribution function for the Poisson
#' binomial distribution. Computational Statistics & Data Analysis, 59, 41-51.
#'
#' @references
#' Volkova, A. Y. (1996). A refinement of the central limit theorem for sums
#' of independent random indicators. Theory of Probability and its Applications 40, 791-794.
#'
#' @references
#' Choi, K.P. and Xia, A. (2002). Approximating the number of successes in independent trials:
#' Binomial versus Poisson. The Annals of Applied Probability, 14(4), 1139-1148.
#'
#' @name PoisBinom
#' @aliases PoisBinom
#' @aliases rpbinom
#' @keywords distribution
#'
#' @export

dpbinom <- function(x, prob, log = FALSE,
                    method = c("MC", "PA", "NA", "BA", "RA", "RNA"),
                    nsim = 10000L) {

  stopifnot(all(prob >= 0 & prob <= 1))
  method <- match.arg(method)

  if (method == "PA") {
    dpois(x, sum(prob), log)
  } else if (method == "NA") {
    dnorm(x, sum(prob), sqrt(sum(prob*(1-prob))), log)
  } else if (method == "BA") {
    dbinom(x, length(prob), mean(prob), log)
  } else {
    if (method == "RNA") {
      p <- .Call('extraDistr_rna_ppbinom', PACKAGE = 'extraDistr', prob)
      for (i in length(p):2)
        p[i] <- p[i]-p[i-1]
    } else {
      if (method == "RA") {
        if (length(prob) >= 20)
          warning("Recursive algorithm is unstable for number of trials of 20 or more.")
        p <- .Call('extraDistr_ra_dpbinom', PACKAGE = 'extraDistr', prob)
      } else {
        p <- .Call('extraDistr_sim_dpbinom', PACKAGE = 'extraDistr', prob, nsim)
      }
    }
    p <- p[x+1]
    if (log) log(p)
    else p
  }
}


#' @rdname PoisBinom
#' @export

ppbinom <- function(x, prob, lower.tail = TRUE, log.p = FALSE,
                    method = c("MC", "PA", "NA", "BA", "RA", "RNA"),
                    nsim = 10000L) {

  stopifnot(all(prob >= 0 & prob <= 1))
  method <- match.arg(method)

  if (method == "PA") {
    ppois(x, sum(prob), lower.tail, log.p)
  } else if (method == "NA") {
    pnorm(x, sum(prob), sqrt(sum(prob*(1-prob))), lower.tail, log.p)
  } else if (method == "BA") {
    pbinom(x, length(prob), mean(prob), lower.tail, log.p)
  } else {
    if (method == "RNA") {
      p <- .Call('extraDistr_rna_ppbinom', PACKAGE = 'extraDistr', prob)
    } else {
      if (method == "RA") {
        if (length(prob) >= 20)
          warning("Recursive algorithm is unstable for number of trials of 20 or more.")
        p <- .Call('extraDistr_ra_dpbinom', PACKAGE = 'extraDistr', prob)
      } else {
        p <- .Call('extraDistr_sim_dpbinom', PACKAGE = 'extraDistr', prob, nsim)
      }
      p <- cumsum(p)
    }
    if (!lower.tail)
      p <- 1-p
    p <- p[x+1]
    if (log.p) log(p)
    else p
  }
}


#' @rdname PoisBinom
#' @export

rpbinom <- function(n, prob) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rpbinom', PACKAGE = 'extraDistr', n, prob)
}

