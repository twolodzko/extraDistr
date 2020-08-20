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


/*
 *  Rayleigh distribution
 *
 *  Values:
 *  x >= 0
 *
 *  Parameters:
 *  sigma > 0
 *
 *  f(x)    = x/sigma^2 * exp(-(x^2 / 2*sigma^2))
 *  F(x)    = 1 - exp(-x^2 / 2*sigma^2)
 *  F^-1(p) = sigma * sqrt(-2 * log(1-p))
 *
 */


inline double logpdf_rayleigh(double x, double sigma,
                              bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(sigma))
    return x+sigma;
#endif
  if (sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x <= 0.0 || !R_FINITE(x))
    return R_NegInf;
  // x/(sigma*sigma) * exp(-(x*x) / (2.0*(sigma*sigma)));
  double lsigsq = 2.0 * log(sigma);
  double lxsq = 2.0 * log(x);
  return log(x) - lsigsq - exp( lxsq - LOG_2F - lsigsq );
}

inline double cdf_rayleigh(double x, double sigma,
                           bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(sigma))
    return x+sigma;
#endif
  if (sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0)
    return 0.0;
  if (!R_FINITE(x))
    return 1.0;
  return 1.0 - exp(-(x*x) / (2.0*(sigma*sigma)));
}

inline double invcdf_rayleigh(double p, double sigma,
                              bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(sigma))
    return p+sigma;
#endif
  if (!VALID_PROB(p) || sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  return sqrt(-2.0*(sigma*sigma) * log(1.0-p));
}

inline double rng_rayleigh(double sigma, bool& throw_warning) {
  if (ISNAN(sigma) || sigma <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  double u = rng_unif();
  return sqrt(-2.0*(sigma*sigma) * log(u));
}


// [[Rcpp::export]]
NumericVector cpp_drayleigh(
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
    p[i] = logpdf_rayleigh(GETV(x, i), GETV(sigma, i),
                           throw_warning);

  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_prayleigh(
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
    p[i] = cdf_rayleigh(GETV(x, i), GETV(sigma, i),
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
NumericVector cpp_qrayleigh(
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
    q[i] = invcdf_rayleigh(GETV(pp, i), GETV(sigma, i),
                           throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rrayleigh(
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
    x[i] = rng_rayleigh(GETV(sigma, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

