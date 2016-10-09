
### 1.8.2

* Corrected and simplified documentation for `*prop` distribution
* Categorical, multinomial, mixture of normals and mixture of Poisson
  distributions are now *less restrictive* about probability parameters
  and accept any non-negative values. Probability parameter vectors are
  normalized to sum to one (i.e. `c(1,1,1)` becomes `c(1,1,1)/3`)
* `NA`'s and `NaN`'s in input now always lead to `NA`'s in output
* Fixed bugs in `rtnorm`, now it properly handles sampling from extreme
  tails

### 1.8.1

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

* Cleaning namespace - now mathematical functions are explicitely
  called from std library
* "using namespace Rcpp" was removed from all the files
* All the numerical values are now explicitely `<double>`'s,
  or casted to `<double>`; or `<int>`'s (for indexing)
* Improvements in algorithms for discrete uniform, categorical,
  mixture of normal, mixture of Poisson distributions
* Improvements in discrete uniform; now it accepts only integer
  valued parameters.
* Major C++ code clean-up

### 1.6.14

* Removed `erf`, `erfc`, `inv_erf` that are not used at this moment and
  caused problems when compiling on Fedora and Solaris

### 1.6.13

* cran-comments.Rd is now in .Rbuildignore

### 1.6.12

* Minor improvements in documentation

### 1.6.11

* Added mixtures of normal and of Poisson distributions

### 1.6.10

* Added truncated Poisson distribution

### 1.6.9

* Minor bug fixes
* Minor improvements in documention
* Additional tests
* Truncated normal returns normal for infinite truncation points

### 1.6.8

* Added half-t, half-normal and half-Cauchy distributions
* Minor changes in documentation
* Changed naming of data-variable from `x` to `q` for CDF's
* Added inverse-CDF for discrete normal distribution, fixed 
  random generation

### 1.6.7

* Improvements in the C++ code

### 1.6.6

* Improvements in documentation
* Added tests
* Multiple minor bug fixes (e.g. functions returning `NaN`
  instead of `0` for `Inf` values)

### 1.6.5

* Added Birnbaum-Saunders (fatigue life) distribution

### 1.6.4

* New algorithm for `rtriang`
* Minor changes in documentation

### 1.6.3

* Added quantile functions for: zero-inflated Poisson, zero-inflated binomial
  zero-inflated negative binomial, inverse gamma, inverse chi-squared
  distributions
* Added Huber density

### 1.6.2

* Minor improvements in the code

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

### 1.5.19

* Ranamed zero-inflated distributions to *`zip` and *`zinb`
* Added zero-inflated binomial *`zib`

### 1.5.18

* Code clean-up
* Added `pzipois` and `pzibinom`

### 1.5.17

* Added `qlgser`; changes in `plgser` and `rlgser`
* Fixed bug in `rlgser`

### 1.5.16

* Examples for most of the functions

### 1.5.15

* Exception handling for `dmvhyper` and `rmvhyper`: values of `x`, `n`,
  and `k` are checked against being non-integers

### 1.5.14

* Added README and NEWS files

