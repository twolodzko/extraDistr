
The package implements several univariate and multivariate distributions
(density, distribution functions, quantile functions, random generation)
that are not available, or are available across multiple other libraries.
The functions are in each case vectorized (there is vectorized version of
multinomial distribution -- unvectorized in base R).

All the functions mimic behavior of the d*, p*, q*, r* functions in R.
Each function has implemented exception handling for improper imput values
and for improper parameter values. For improper input values density or
probability mass functions returns 0's. For improper parameter values NaN's
are returned with warning. For discrete distributions, 0's are returned by
probability mass functions for non-discrete values with warning.

Checks using devtools check() and build_win() yield no errors or warnings
besides the following note appears on linux:

* checking installed package size ... NOTE
  installed size is  8.9Mb
  sub-directories of 1Mb or more:
    libs   8.6Mb
    
It is similar to the one described recently [1]. This package depends
only on Rcpp and some standard linux libraries (stl) and does not use
any redundant C++ libraries, so it's size is adequate to it's pourpose.

Behavior of the functions was tested using testthat package (tests/)
and additional checks are available through examples. Most of the random
generation is done by inverse transform, so valid quantile function leads
to valid random generation.

In several cases log's are used to prevent underflow. When summing
discrete probabilities they are multiplied by large constant
(P_NORM_CONST in const.h) for the same reason. When checking if
sum(prob) == 1 (e.g. multinomial distribution) values are compared using
predefined tolerance (MIN_DIFF_EPS in const.h) rather then checking for
equality.


[1]: https://borishejblum.wordpress.com/2016/04/27/cran-check-note-sub-directories-of-1mb-or-more-libs/