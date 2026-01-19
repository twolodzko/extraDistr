
# ⚠️ This project has changed ownership

After several years of maintaining it, I decided to retire and pass over the maintainance of it to the community.
You can find the [contact to the new maintainer(s) on CRAN](https://cran.r-project.org/web/packages/extraDistr/index.html).

---

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
