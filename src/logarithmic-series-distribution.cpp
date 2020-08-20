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

/*
*  Logarithmic Series distribution
*
*  Values:
*  x
*
*  Parameters:
*  0 < theta < 1
*
*  f(x) = (-1/log(1-theta)*theta^x) / x
*  F(x) = -1/log(1-theta) * sum((theta^x)/x)
*
*/


inline double logpdf_lgser(double x, double theta, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(theta))
    return x+theta;
#endif
  if (theta <= 0.0 || theta >= 1.0) {
    throw_warning = true;
    return NAN;
  }
  if (!isInteger(x) || x < 1.0)
    return R_NegInf;
  // a = -1.0/log(1.0 - theta);
  double a = -1.0/log1p(-theta);
  // a * pow(theta, x) / x;
  return log(a) + (log(theta) * x) - log(x);
}

inline double cdf_lgser(double x, double theta, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(theta))
    return x+theta;
#endif
  if (theta <= 0.0 || theta >= 1.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 1.0)
    return 0.0;
  if (!R_FINITE(x))
    return 1.0;
  if (is_large_int(x)) {
    Rcpp::warning("NAs introduced by coercion to integer range");
    return NA_REAL;
  }
  
  double a = -1.0/log1p(-theta);
  double b = 0.0;
  double dk;
  int ix = to_pos_int(x);
  
  for (int k = 1; k <= ix; k++) {
    dk = to_dbl(k);
    b += pow(theta, dk) / dk;
  }
  
  return a * b;
}

inline double invcdf_lgser(double p, double theta, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(theta))
    return p+theta;
#endif
  if (theta <= 0.0 || theta >= 1.0 || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  if (p == 0.0)
    return 1.0;
  if (p == 1.0)
    return R_PosInf;
  
  double pk = -theta/log(1.0 - theta);
  double k = 1.0;
  
  while (p > pk) {
    p -= pk;
    pk *= theta * k/(k+1.0);
    k += 1.0;
  }
  
  return k;
}

inline double rng_lgser(double theta, bool& throw_warning) {
  if (ISNAN(theta) || theta <= 0.0 || theta >= 1.0) {
    throw_warning = true;
    return NA_REAL;
  }

  double u = rng_unif();
  double pk = -theta/log(1.0 - theta);
  double k = 1.0;
  
  while (u > pk) {
    u -= pk;
    pk *= theta * k/(k+1.0);
    k += 1.0;
  }
  
  return k;
}


// [[Rcpp::export]]
NumericVector cpp_dlgser(
    const NumericVector& x,
    const NumericVector& theta,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), theta.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    theta.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_lgser(GETV(x, i), GETV(theta, i),
                        throw_warning);
 
 if (!log_prob)
   p = Rcpp::exp(p);
 
 if (throw_warning)
   Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_plgser(
    const NumericVector& x,
    const NumericVector& theta,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), theta.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    theta.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_lgser(GETV(x, i), GETV(theta, i),
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
NumericVector cpp_qlgser(
    const NumericVector& p,
    const NumericVector& theta,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), theta.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    p.length(),
    theta.length()
  });
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_lgser(GETV(pp, i), GETV(theta, i),
                        throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rlgser(
    const int& n,
    const NumericVector& theta
  ) {
  
  if (theta.length() < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }

  NumericVector x(n);
  
  bool throw_warning = false;

  for (int i = 0; i < n; i++)
    x[i] = rng_lgser(GETV(theta, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

