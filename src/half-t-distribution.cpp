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
 * 
 * x >= 0
 * 
 * Parameters:
 * nu > 0
 * sigma > 0
 * 
 * with nu = 1   returns half-Cauchy
 * with nu = Inf returns half-normal
 * 
 */

inline double pdf_ht(double x, double nu, double sigma,
                     bool& throw_warning) {
  if (ISNAN(x) || ISNAN(nu) || ISNAN(sigma))
    return x+nu+sigma;
  if (sigma <= 0.0 || nu <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  return 2.0 * R::dt(x/sigma, nu, false)/sigma;
}

inline double cdf_ht(double x, double nu, double sigma,
                     bool& throw_warning) {
  if (ISNAN(x) || ISNAN(nu) || ISNAN(sigma))
    return x+nu+sigma;
  if (sigma <= 0.0 || nu <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  return 2.0 * R::pt(x/sigma, nu, true, false) - 1.0;
}

inline double invcdf_ht(double p, double nu, double sigma,
                        bool& throw_warning) {
  if (ISNAN(p) || ISNAN(nu) || ISNAN(sigma))
    return p+nu+sigma;
  if (sigma <= 0.0 || nu <= 0.0 || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  return R::qt((p+1.0)/2.0, nu, true, false) * sigma;
}

inline double rng_ht(double nu, double sigma, bool& throw_warning) {
  if (ISNAN(nu) || ISNAN(sigma) || sigma <= 0.0 || nu <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  return abs(R::rt(nu) * sigma);
}


// [[Rcpp::export]]
NumericVector cpp_dht(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& sigma,
    const bool& log_prob = false
  ) {
  
  int Nmax = std::max({
    x.length(),
    nu.length(),
    sigma.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++) 
    p[i] = pdf_ht(GETV(x, i), GETV(nu, i),
                  GETV(sigma, i), throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pht(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  int Nmax = std::max({
    x.length(),
    nu.length(),
    sigma.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_ht(GETV(x, i), GETV(nu, i),
                  GETV(sigma, i), throw_warning);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qht(
    const NumericVector& p,
    const NumericVector& nu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  int Nmax = std::max({
    p.length(),
    nu.length(),
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
    q[i] = invcdf_ht(GETV(pp, i), GETV(nu, i),
                     GETV(sigma, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rht(
    const int& n,
    const NumericVector& nu,
    const NumericVector& sigma
  ) {
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_ht(GETV(nu, i), GETV(sigma, i),
                  throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

