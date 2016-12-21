
# extraDistr

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/twolodzko/extraDistr.svg?branch=master)](https://travis-ci.org/twolodzko/extraDistr)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/extraDistr)](http://cran.r-project.org/package=extraDistr)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/extraDistr?color=orange)](http://cranlogs.r-pkg.org/badges/grand-total/extraDistr)

Density, distribution function, quantile function and random
generation for a number of univariate and multivariate distributions.

This package follows naming convention that is consistent with base R,
where density (or probability mass) functions, distribution functions,
quantile functions and random generation functions names are followed by
`d`\*, `p`\*, `q`\*, and `r`\* prefixes.

Behaviour of the functions mimics the base R, where for
invalid parameters `NaN`'s are returned, while
for values beyond function support 0's are returned
(e.g. for non-integers in discrete distributions, or for
negative values in functions with non-negative support).

All the functions vectorized and coded in C++11 using Rcpp.
