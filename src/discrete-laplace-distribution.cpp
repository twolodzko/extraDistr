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


inline double pmf_dlaplace(double x, double p, double mu,
                           bool& throw_warning) {
  if (ISNAN(x) || ISNAN(p) || ISNAN(mu))
    return x+p+mu;
  if (p <= 0.0 || p >= 1.0) {
    throw_warning = true;
    return NAN;
  }
  if (!isInteger(x))
    return 0.0;
  return (1.0-p)/(1.0+p) * pow(p, abs(x-mu));
} 

inline double cdf_dlaplace(double x, double p, double mu,
                           bool& throw_warning) {
  if (ISNAN(x) || ISNAN(p) || ISNAN(mu))
    return x+p+mu;
  if (p <= 0.0 || p >= 1.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return pow(p, -floor(x-mu))/(1.0+p);
  else
    return 1.0 - (pow(p, floor(x-mu)+1.0)/(1.0+p));
} 

inline double rng_dlaplace(double p, double mu,
                           bool& throw_warning) {
  if (ISNAN(p) || ISNAN(mu) || p <= 0.0 || p >= 1.0) {
    throw_warning = true;
    return NA_REAL;
  }
  double q, u, v;
  q = 1.0 - p;
  u = R::rgeom(q); 
  v = R::rgeom(q); 
  return u-v + mu;
} 


// [[Rcpp::export]]
NumericVector cpp_ddlaplace(
    const NumericVector& x,
    const NumericVector& scale,
    const NumericVector& location,
    const bool& log_prob = false
) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(scale.length());
  dims.push_back(location.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pmf_dlaplace(GETV(x, i), GETV(scale, i),
                        GETV(location, i), throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pdlaplace(
    const NumericVector& x,
    const NumericVector& scale,
    const NumericVector& location,
    const bool& lower_tail = true,
    const bool& log_prob = false
) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(scale.length());
  dims.push_back(location.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_dlaplace(GETV(x, i), GETV(scale, i),
                        GETV(location, i), throw_warning);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rdlaplace(
    const int& n,
    const NumericVector& scale,
    const NumericVector& location
) {
  
  std::vector<int> dims;
  dims.push_back(scale.length());
  dims.push_back(location.length());
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_dlaplace(GETV(scale, i), GETV(location, i),
                        throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}
