
#' Additional univariate and multivariate distributions
#'
#' Density, distribution function, quantile function and random
#' generation for a number of univariate and multivariate distributions.
#' 
#' @details
#' 
#' This package follows naming convention that is consistent with base R,
#' where density (or probability mass) functions, distribution functions,
#' quantile functions and random generation functions names are followed by
#' \code{d}*, \code{p}*, \code{q}*, and \code{r}* prefixes.
#' 
#' Behaviour of the functions is consistent with base R, where for
#' not valid parameters values \code{NaN}'s are returned, while
#' for values beyond function support \code{0}'s are returned
#' (e.g. for non-integers in discrete distributions, or for
#' negative values in functions with non-negative support).
#' 
#' All the functions vectorized and coded in C++ using \pkg{Rcpp}.
#' 
#' @docType package
#' @name extraDistr
#' 
#' @useDynLib extraDistr, .registration = TRUE 
#' @importFrom Rcpp sourceCpp
#' 
#' @importFrom stats pgamma qgamma rgamma
#' @importFrom stats pnorm qnorm rnorm
#' @importFrom stats rpois
#' @importFrom stats qchisq rchisq
NULL
