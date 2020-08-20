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
* Zero-inflated Poisson distribution
* 
* Parameters:
* lambda > 0
* 0 <= pi <= 1
* 
* Values:
* x >= 0
*
*/

inline double pdf_zip(double x, double lambda, double pi,
                      bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(pi))
    return x+lambda+pi;
#endif
  if (lambda <= 0.0 || !VALID_PROB(pi)) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0 || !isInteger(x) || !R_FINITE(x))
    return 0.0;
  if (x == 0.0) {
    // pi + (1.0-pi) * exp(-lambda);
    return pi + exp( log1p(-pi) - lambda );
  } else {
    // (1.0-pi) * R::dpois(x, lambda, false);
    return exp( log1p(-pi) + R::dpois(x, lambda, true) );
  }
}

inline double cdf_zip(double x, double lambda, double pi,
                      bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(pi))
    return x+lambda+pi;
#endif
  if (lambda <= 0.0 || !VALID_PROB(pi)) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  if (!R_FINITE(x))
    return 1.0;
  // pi + (1.0-pi) * R::ppois(x, lambda, true, false);
  return pi + exp(log1p(-pi) + R::ppois(x, lambda, true, true));
}

inline double invcdf_zip(double p, double lambda, double pi,
                         bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(lambda) || ISNAN(pi))
    return p+lambda+pi;
#endif
  if (lambda <= 0.0 || !VALID_PROB(pi) || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  if (p < pi)
    return 0.0;
  else
    return R::qpois((p - pi) / (1.0-pi), lambda, true, false);
}

inline double rng_zip(double lambda, double pi, bool& throw_warning) {
  if (ISNAN(lambda) || ISNAN(pi) ||
      lambda <= 0.0 || !VALID_PROB(pi)) {
    throw_warning = true;
    return NA_REAL;
  }
  double u = rng_unif();
  if (u < pi)
    return 0.0;
  else
    return R::rpois(lambda);
}


// [[Rcpp::export]]
NumericVector cpp_dzip(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& pi,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), lambda.length(), pi.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    lambda.length(),
    pi.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_zip(GETV(x, i), GETV(lambda, i),
                   GETV(pi, i), throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pzip(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& pi,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), lambda.length(), pi.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    lambda.length(),
    pi.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_zip(GETV(x, i), GETV(lambda, i),
                   GETV(pi, i), throw_warning);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qzip(
    const NumericVector& p,
    const NumericVector& lambda,
    const NumericVector& pi,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), lambda.length(), pi.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    p.length(),
    lambda.length(),
    pi.length()
  });
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_zip(GETV(pp, i), GETV(lambda, i),
                      GETV(pi, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rzip(
    const int& n,
    const NumericVector& lambda,
    const NumericVector& pi
  ) {
  
  if (std::min({lambda.length(), pi.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_zip(GETV(lambda, i), GETV(pi, i),
                   throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return x;
}

