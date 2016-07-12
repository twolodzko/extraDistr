

#' Discrete Laplace distribution
#' 
#' Probability mass, distribution function and random generation
#' for the discrete Laplace distribution parametrized by location and scale.
#' 
#' @param x,q             vector of quantiles.
#' @param n	              number of observations. If \code{length(n) > 1},
#'                        the length is taken to be the number required.
#' @param scale           scale parameter; \code{0 < scale < 1}.
#' @param location        location parameter.
#' @param log,log.p	      logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	    logical; if TRUE (default), probabilities are \eqn{P[X \le x]}
#'                        otherwise, \eqn{P[X > x]}.
#' 
#' @details 
#' 
#' If \eqn{U \sim \mathrm{Geometric}(1-p)}{U ~ Geometric(1-p)} and
#' \eqn{V \sim \mathrm{Geometric}(1-p)}{V ~ Geometric(1-p)},
#' then \eqn{U-V \sim \mathrm{DiscreteLaplace}(p)}{U-V ~ DiscreteLaplace(p)},
#' where geometric distribution is related to discrete Laplace distribution
#' in similar way as exponential distribution is related to Laplace distribution.
#' 
#' Probability mass function
#' 
#' \deqn{
#' f(x) = \frac{1-p}{1+p} p^{|x-\mu|}
#' }{
#' f(x) = (1-p)/(1+p) * p^(|x-\mu|)
#' }
#' 
#' Cumulative distribution function
#' 
#' \deqn{
#' F(x) = \left\{\begin{array}{ll}
#' \frac{p^{-|x-\mu|}}{1+p} & x < 0 \\
#' 1 - \frac{p^{|x-\mu|+1}}{1+p} & x \ge 0
#' \end{array}\right.
#' }{
#' F(x) = [if x < 0:] p^-floor(x-\mu))/(1+p) [else:] 1-(p^(floor(x-\mu)+1))/(1+p)
#' }
#' 
#' @references 
#' Inusah, S., & Kozubowski, T.J. (2006). A discrete analogue of the Laplace distribution.
#' Journal of statistical planning and inference, 136(3), 1090-1102.
#' 
#' @references 
#' Kotz, S., Kozubowski, T., & Podgorski, K. (2012).
#' The Laplace distribution and generalizations: a revisit with applications
#' to communications, economics, engineering, and finance.
#' Springer Science & Business Media.
#' 
#' @examples 
#' 
#' p <- 0.45
#' x <- rdlaplace(1e5, p)
#' xx <- seq(-200, 200, by = 1)
#' plot(prop.table(table(x)))
#' lines(xx, ddlaplace(xx, p), col = "red")
#' plot(hist(pdlaplace(x, p)))
#' plot(ecdf(x))
#' lines(xx, pdlaplace(xx, p), col = "red")
#' 
#' @name DiscreteLaplace
#' @aliases DiscreteLaplace
#' @aliases ddlaplace
#' @keywords distribution
#' 
#' @export

ddlaplace <- function(x, scale, location = 0, log = FALSE) {
  .Call('extraDistr_cpp_ddlaplace', PACKAGE = 'extraDistr', x, scale, location, log)
}


#' @rdname DiscreteLaplace
#' @export

pdlaplace <- function(q, scale, location = 0, lower.tail = TRUE, log.p = FALSE) {
  .Call('extraDistr_cpp_pdlaplace', PACKAGE = 'extraDistr', q, scale, location, lower.tail, log.p)
}


#' @rdname DiscreteLaplace
#' @export

rdlaplace <- function(n, scale, location = 0) {
  if (length(n) > 1) n <- length(n)
  .Call('extraDistr_cpp_rdlaplace', PACKAGE = 'extraDistr', n, scale, location)
}



# ==============================================
# 
# ddlaplace <- function(x, p, mu = 0) (1-p)/(1+p)*p^abs(x-mu)
# pdlaplace <- function(q, p, mu = 0) ifelse(q<0, (p^-floor(q-mu))/(1+p), 1-(p^(floor(q-mu)+1))/(1+p))
# rdlaplace <- function(n, p, mu = 0) rgeom(n, 1-p) - rgeom(n, 1-p) + mu
# 
# qdlaplace <- function(pp, p, mu = 0) {
#   qq <- (pp - 0.5) * 2
#   ifelse(qq>0, qgeom(qq, 1-p), -qgeom(abs(qq), 1-p))
# }
# 
# p <- 0.46
# x <- rdlaplace(1e5, p)
# pp <- seq(0, 1, by = 0.0001)
# plot(ecdf(x))
# lines(qdlaplace(pp, p), pp, col = "orange")



# # ==============================================
# 
# 
# qdlaplace <- function(pp, p, mu = 0) {
#   qq <- (pp - 0.5) * 2
#   ifelse(qq>0, qgeom(qq, 1-p), -qgeom(abs(qq), 1-p)) 
# }
# 
# x <- rdlaplace(1e5, p)
# pp <- seq(0, 1, by = 0.0001)
# plot(ecdf(x))
# lines(qdlaplace(pp, p), pp, col = "orange")
# 
# # 
# # grint <- function(x) ifelse(x>0, floor(x), ceiling(x))
# # 
# # xx <- seq(-30, 30, by = 1)
# # 
# # plot(x, f(x, p = 0.5), type = "b", ylim = c(0, 0.9))
# # lines(x, f(x, p = 0.3), type = "b", col = "red")
# # lines(x, f(x, p = 0.1), type = "b", col = "orange")
# # lines(x, f(x, p = 0.6), type = "b", col = "blue")
# # lines(x, f(x, p = 0.8), type = "b", col = "violet")
# # 
# # xx <- seq(-200, 200, by = 1)
# # 
# # p <- 0.45
# # # x1 <- rgeom(1e5, p = p)
# # # x2 <- rgeom(1e5, p = p)
# # # 
# # # y <- x1-x2
# # x <- rdlaplace(1e5, p)
# # plot(prop.table(table(x)))
# # lines(xx, ddlaplace(xx, p), col = "red")
# # 
# # plot(ecdf(x))
# # lines(xx, pdlaplace(xx, p), col = "red")
# # 
# # plot(hist(pdlaplace(x, p)))
