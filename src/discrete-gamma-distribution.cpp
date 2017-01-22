#include <Rcpp.h>
#include "shared.h"

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;


/*
* Discrete normal distribution
* 
* Values:
* x
* 
* Parameters
* mu
* sigma > 0
*  
*/


inline double pmf_dgamma(double x, double shape, double scale,
                         bool& throw_warning) {
  if (ISNAN(x) || ISNAN(shape) || ISNAN(scale))
    return x+shape+scale;
  if (shape <= 0.0 || scale <= 0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0 || !isInteger(x))
    return 0.0;
  return R::pgamma(x+1.0, shape, scale, true, false) -
    R::pgamma(x, shape, scale, true, false);
}


// [[Rcpp::export]]
NumericVector cpp_ddgamma(
    const NumericVector& x,
    const NumericVector& shape,
    const NumericVector& scale,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(shape.length());
  dims.push_back(scale.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pmf_dgamma(GETV(x, i), GETV(shape, i),
                      GETV(scale, i), throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}

