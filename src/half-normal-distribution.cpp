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


inline double logpdf_hnorm(double x, double sigma,
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
  return LOG_2F + R::dnorm(x, 0.0, sigma, true);
}

inline double cdf_hnorm(double x, double sigma, bool& throw_warning) {
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
  return 2.0 * R::pnorm(x, 0.0, sigma, true, false) - 1.0;
}

inline double invcdf_hnorm(double p, double sigma, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(sigma))
    return p+sigma;
#endif
  if (sigma <= 0.0 || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  return R::qnorm((p+1.0)/2.0, 0.0, sigma, true, false);
}

inline double rng_hnorm(double sigma, bool& throw_warning) {
  if (ISNAN(sigma) || sigma <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  return abs(R::norm_rand()) * sigma;
}


// [[Rcpp::export]]
NumericVector cpp_dhnorm(
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
    p[i] = logpdf_hnorm(GETV(x, i), GETV(sigma, i),
                        throw_warning);
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_phnorm(
    const NumericVector& x,
    const NumericVector& sigma,
    const bool& lower_tail = true,
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
    p[i] = cdf_hnorm(GETV(x, i), GETV(sigma, i),
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
NumericVector cpp_qhnorm(
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
    q[i] = invcdf_hnorm(GETV(pp, i), GETV(sigma, i),
                        throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rhnorm(
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
    x[i] = rng_hnorm(GETV(sigma, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

