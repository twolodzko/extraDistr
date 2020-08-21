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
using std::tan;
using std::atan;


inline double logpdf_hcauchy(double x, double sigma,
                             bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(sigma))
    return x+sigma;
#endif
  if (sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return R_NegInf;
  // 2.0/(M_PI*(1.0 + pow(x/sigma, 2.0)))/sigma;
  return LOG_2F - log(M_PI) - log1p(exp( (log(x)-log(sigma)) * 2.0 )) - log(sigma);
}

inline double cdf_hcauchy(double x, double sigma, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(sigma))
    return x+sigma;
#endif
  if (sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  return 2.0/M_PI * atan(x/sigma);
}

inline double invcdf_hcauchy(double p, double sigma,
                             bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(sigma))
    return p+sigma;
#endif
  if (sigma <= 0.0 || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  return sigma * tan((M_PI*p)/2.0);
}

inline double rng_hcauchy(double sigma, bool& throw_warning) {
  if (ISNAN(sigma) || sigma <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  return abs(R::rcauchy(0.0, sigma));
}


// [[Rcpp::export]]
NumericVector cpp_dhcauchy(
    const NumericVector& x,
    const NumericVector& sigma,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    sigma.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_hcauchy(GETV(x, i), GETV(sigma, i),
                          throw_warning);
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_phcauchy(
    const NumericVector& x,
    const NumericVector& sigma,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  if (std::min({x.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    sigma.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_hcauchy(GETV(x, i), GETV(sigma, i),
                       throw_warning);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qhcauchy(
    const NumericVector& p,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    p.length(),
    sigma.length()
  });
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_hcauchy(GETV(pp, i), GETV(sigma, i),
                          throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rhcauchy(
    const int& n,
    const NumericVector& sigma
  ) {
  
  if (sigma.length() < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_hcauchy(GETV(sigma, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

