### 1.10.0

* Fixed bug in `rgpd` which produced negative samples.
* The RcppExports.cpp was fixed by reformatting it thanks to Dirk Eddelbuettel.

### 1.9.1

* Generated header file, `inst/include/extraDistr.h`, to make C++ code callable 
  from C++ code in other R packages (#23, contribution from Jon Fintzi, @fintzij).

### 1.8.11

* Fixed bugs in Frechet distribution implementations (`pfrechet`).

### 1.8.10

* Now zero-inflated negative-binomial (#12) and beta negative-binomial distributions (#14)
  do not have integer-only constraint on `size` parameter.
* Fixed mistakes in the zero-inflated Poisson distribution documentation (#15) and
  negative hypergeometric distribution documentation.
* Minor changes towards future Rcpp STRICT_R_HEADERS compatibility.
* Fixed bug in vectorization code for `pbbinom` (#16). Additionally, this should make 
  `pbbinom` and `pbnbinom` faster when working with large vectors.

### 1.8.9

* Fixed bug in `dinvgamma` function
* Updated the DESCRIPTION file by mentioning the packages used in
  optional unit tests for CRAN compatibility

### 1.8.8

* Deprecated functions: `dnst`, `pnst`, `qnst` and `rnst` were removed
* Fixed typos in documentation (thanks to #8 by philchalmers)
* Converted a number of pdf and cdf functions to more numerically stable
  versions using logs (see #9 and #10)
* Fixed bug in `ppower`: with `lower.tail = FALSE` it returned wrong values
* Fixed bug in `dgpd` and `pgpd`: they assumed slightly wrong support 
* Improved the `rgev` and `rgpd`, now they give better random values since 
  relying on exponential distribution random generator if possible
* Documentation was improved and corrected in several places
* Power distribution functions now check if `alpha > 0` and `beta > 0`

### 1.8.7

* Fixed bug in `pinvgamma` (`lower.tail` and `log.p` didn't work)
* Fixed underflow issues with `rmnom` and `rdirmnom` (#7)
* The `*nst` functions are now deprecated and renamed to more informative
  abbrevation `*lst`

### 1.8.6

* Now consistently with base R only the first elements of the logical
  arguments are used (thanks to #5)
* Fixed bug in `rtnorm` (sampling from lower bound returned incorrect
  values)
* When computation becomes slow, now `pbnbinom`, `pbbinom`, `pgpois`
  functions are easier to brake
* Automatically registering native routines via Rcpp
* Fixed bug in `pinvgamma` (it returned non-zero probabilities for
  q < 0)
  
### 1.8.5

* Now `rmnom` and `rdirmnom` (issue #3) do not return `NaN`'s due to
  underflow issues
* Fixed bug (issue #4) that resulted in hanging R if zero-length vectors
  were provided as input

### 1.8.4

* Random generation from categorical distribution parametrized by
  log-probabilities `rcatlp`
* Re-parametrized beta distribution is now more flexible since
  "prior" parameter was introduced (see documentation); this change
  is connected to switching to different prior equal to zero
  (instead of one as previously). Such parametrization is more
  commonly used in the literature. This change is documented in
  the `?PropBeta` documentation
* Improvements in documentation and examples
* Fixed bugs in `pbbinom`, `pbnbinom`, `pgpois`, `*nhyper` that
  prevented compiling on RedHat Linux (#2)

### 1.8.3

* Switched to C++11, underlying code was simplified and improved
* Using memoization techniques for `pbbinom`, `pbnbinom`, `pgpois` and
  negative hypergeometric distribution that lead to major speed improvements
* Improved and simplified code for beta-binomial, beta negative-binomial,
  and gamma-Poisson cdf; now recursive algorithm employing memoization techniques
  is used what noticeably improves performance
* Discrete gamma, shifted Gompertz (d,p,r), negative hypergeometric and
  truncated binomial distributions (d,p,q,r) were added
* Now `lower.tail` and `log.p` options for `pbetapr` work properly
* Simplified code for multivariate hypergeometric, multinomial,
  Dirichlet-multinomial and categorical distributions
* Truncated poisson distribution is now parametrized by lower and upper
  truncation points
* Fixed bugs in `dbvpois` (it returned underestimated probabilities),
  `dslash` (there was discontinuity at x=0), `pcat` (randomly it broke
  if x was greater then the upper limit), and `pdnorm`.
* Random generation functions throw warnings and produce `NA`'s on `NA`'s in
  parameters or incorrect parameters - as in base R
* Order of parameters in discrete Laplace distribution was changed to
  location and scale (vs scale and location) for consistency with other
  distributions (e.g. continous Laplace)
* Improved exception handling


### 1.8.1-2

* Corrected and simplified documentation for `*prop` distribution
* Categorical, multinomial, mixture of normals and mixture of Poisson
  distributions are now *less restrictive* about probability parameters
  and accept any non-negative values. Probability parameter vectors are
  normalized to sum to one (i.e. `c(1,1,1)` becomes `c(1,1,1)/3`)
* `NA`'s and `NaN`'s in input now always lead to `NA`'s in output
* Fixed bugs in `rtnorm`, now it properly handles sampling from extreme
  tails
* Tests for multivariate distribution now done
  with tolerance `1e-2`

### 1.8.0

* Bug fixes in `qtnorm` - now it works properly for
  non-standard truncated normal
* Bug fixes for `rmnom` and `rdirmnom` - now they correctly
  use the `prob` parameter values
* Minor improvements in `dmvhyper` and `dbvpois`
* Cleaned-up the documentation
* Added discrete Laplace distribution
* Faster RNG generator for Laplace distribution
* Changes to using lower level RNG functions
  (`unif_rand`, `norm_rand`) when possible
* More tests

### 1.7.2

* Documentation clean-up
* Minor improvements and simplifications in C++ code
* For `min == max` discrete uniform distribution behaves now
  as degenerate distribution

### 1.7.1

* Added Dirichlet-multinomial and beta prime distributions
* Code clean-up for categorical distribution
* Improvements in documentation and examples (e.g. bivariate Poisson,
  bivariate normal)
* Other minor improvements in the code

### 1.7.0

* Cleaning namespace - now mathematical functions are explicitly
  called from std library
* "using namespace Rcpp" was removed from all the files
* All the numerical values are now explicitly `<double>`'s,
  or casted to `<double>`; or `<int>`'s (for indexing)
* Improvements in algorithms for discrete uniform, categorical,
  mixture of normal, mixture of Poisson distributions
* Improvements in discrete uniform; now it accepts only integer
  valued parameters
* Major C++ code clean-up

### 1.6.10-14

* Removed `erf`, `erfc`, `inv_erf` that are not used at this moment and
  caused problems when compiling on Fedora and Solaris
* Minor improvements in documentation
* Added mixtures of normal and of Poisson distributions
* Added truncated Poisson distribution

### 1.6.7-9

* Minor changes and improvements in documention
* Truncated normal returns normal for infinite truncation points
* Added half-t, half-normal and half-Cauchy distributions
* Changed naming of data-variable from `x` to `q` for CDF's
* Added inverse-CDF for discrete normal distribution, fixed 
  random generation
* Added tests
* Multiple minor bug fixes (e.g. functions returning `NaN`
  instead of `0` for `Inf` values)
* Added Birnbaum-Saunders (fatigue life) distribution
* New algorithm for `rtriang`
* Added quantile functions for: zero-inflated Poisson, zero-inflated binomial
  zero-inflated negative binomial, inverse gamma, inverse chi-squared
  distributions
* Added Huber density

### 1.6.1

* Exception handling in discrete uniform distribution functions
* Bug fixes in `rdunif`
* Clean-up in documentation and examples
* Now `dbvnorm` and `dbvpois` work with matrixes
* Renaming of parameters in bivariate normal to be consistent with
  base R normal distribution

### 1.6.0

* Warning messages for non-integer values in discrete distributions
* Warning when returning NANs
* Code clean-up in `pbbinom` and `pbnbinom`
* Clean-up in documentation

### 1.5.15-19

* Ranamed zero-inflated distributions to `*zip` and `*zinb`
* Added zero-inflated binomial `*zib`
* Code clean-up
* Added `pzipois` and `pzibinom`
* Added `qlgser`; changes in `plgser` and `rlgser`
* Fixed bug in `rlgser`
* Examples for most of the functions
* Exception handling for `dmvhyper` and `rmvhyper`: values of `x`, `n`,
  and `k` are checked against being non-integers

