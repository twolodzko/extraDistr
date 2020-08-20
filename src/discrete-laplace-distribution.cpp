#include <Rcpp.h>
#include "shared.h"
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;

using std::log1p;


inline double logpmf_dlaplace(double x, double p, double mu,
                              bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(p) || ISNAN(mu))
    return x+p+mu;
#endif
  if (p <= 0.0 || p >= 1.0) {
    throw_warning = true;
    return NAN;
  }
  if (!isInteger(x))
    return R_NegInf;
  // (1.0-p)/(1.0+p) * pow(p, abs(x-mu));
  return log1p(-p) - log1p(p) + log(p) * abs(x-mu);
} 

inline double cdf_dlaplace(double x, double p, double mu,
                           bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(p) || ISNAN(mu))
    return x+p+mu;
#endif
  if (p <= 0.0 || p >= 1.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0) {
    // pow(p, -floor(x-mu))/(1.0+p);
    return exp( (log(p) * -floor(x-mu)) - log1p(p) );
  } else {
    // 1.0 - (pow(p, floor(x-mu)+1.0)/(1.0+p))
    return 1.0 - exp( log(p) * (floor(x-mu)+1.0) - log1p(p) );
  }
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
    const NumericVector& location,
    const NumericVector& scale,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), location.length(), scale.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    scale.length(),
    location.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = logpmf_dlaplace(GETV(x, i), GETV(scale, i),
                           GETV(location, i), throw_warning);
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pdlaplace(
    const NumericVector& x,
    const NumericVector& location,
    const NumericVector& scale,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), location.length(), scale.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    scale.length(),
    location.length()
  });
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
    const NumericVector& location,
    const NumericVector& scale
  ) {
  
  if (std::min({location.length(), scale.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_dlaplace(GETV(scale, i), GETV(location, i),
                        throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}
