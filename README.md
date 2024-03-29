
# extraDistr

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/extraDistr)](https://CRAN.R-project.org/package=extraDistr)
[![GitHub Actions CI](https://github.com/twolodzko/extraDistr/workflows/CI/badge.svg)](https://github.com/twolodzko/extraDistr/actions?query=workflow%3ACI)
[![Coverage Status](https://img.shields.io/codecov/c/github/twolodzko/extraDistr/master.svg)](https://app.codecov.io/github/twolodzko/extraDistr?branch=master)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/extraDistr)](https://cran.rstudio.com/web/packages/extraDistr/index.html)


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

All the functions vectorized and coded in C++11 using [Rcpp](https://www.rcpp.org/).
