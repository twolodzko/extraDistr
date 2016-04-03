This is a **development version** of R package implementing a number of
univariate and multivariate probability distributions. It still needs multiple
minor changes and bug fixes - use it at your own risk!

Density, distribution function, quantile function and random
generation for a number of univariate and multivariate distributions.

This package follows naming convention that is consistent with base R,
where density (or probability mass) functions, distribution functions,
quantile functions and random generation functions names are followed by
`d`\*, `p`\*, `q`\*, and `r`\* prefixes.

Behaviour of the functions is consistent with base R, where for
not valid parameters values `NaN`'s are returned, while
for values beyond function support 0's are returned
(e.g. for non-integers in discrete distributions, or for
negative values in functions with non-negative support).

All the functions vectorized and coded in C++ using Rcpp.
