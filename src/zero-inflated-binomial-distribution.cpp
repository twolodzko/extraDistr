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

inline double pdf_zib(double x, double n, double p,
                      double pi, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(n) || ISNAN(p) || ISNAN(pi))
    return x+n+p+pi;
#endif
  if (!VALID_PROB(p) || n < 0.0 || !VALID_PROB(pi) ||
      !isInteger(n, false)) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0 || !isInteger(x) || !R_FINITE(x))
    return 0.0;
  if (x == 0.0) {
    // pi + (1.0-pi) * pow(1.0-p, n);
    return pi + exp( log1p(-pi) + log1p(-p) * n );
  } else {
    // (1.0-pi) * R::dbinom(x, n, p, false);
    return exp( log1p(-pi) + R::dbinom(x, n, p, true) );
  }
}

inline double cdf_zib(double x, double n, double p,
                      double pi, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(n) || ISNAN(p) || ISNAN(pi))
    return x+n+p+pi;
#endif
  if (!VALID_PROB(p) || n < 0.0 || !VALID_PROB(pi) ||
      !isInteger(n, false)) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  if (!R_FINITE(x))
    return 1.0;
  // pi + (1.0-pi) * R::pbinom(x, n, p, true, false);
  return pi + exp( log1p(-pi) + R::pbinom(x, n, p, true, true) );
}

inline double invcdf_zib(double pp, double n, double p,
                         double pi, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(pp) || ISNAN(n) || ISNAN(p) || ISNAN(pi))
    return pp+n+p+pi;
#endif
  if (!VALID_PROB(p) || n < 0.0 || !VALID_PROB(pi) ||
      !isInteger(n, false) || !VALID_PROB(pp)) {
      throw_warning = true;
    return NAN;
  }
  if (pp < pi)
    return 0.0;
  else
    return R::qbinom((pp - pi) / (1.0-pi), n, p, true, false);
}

inline double rng_zib(double n, double p, double pi,
                      bool& throw_warning) {
  if (ISNAN(n) || ISNAN(p) || ISNAN(pi) || !VALID_PROB(p) ||
      n < 0.0 || !VALID_PROB(pi) || !isInteger(n, false)) {
    throw_warning = true;
    return NA_REAL;
  }
  double u = rng_unif();
  if (u < pi)
    return 0.0;
  else
    return R::rbinom(n, p);
}


// [[Rcpp::export]]
NumericVector cpp_dzib(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), size.length(),
                prob.length(), pi.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    size.length(),
    prob.length(),
    pi.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_zib(GETV(x, i), GETV(size, i),
                   GETV(prob, i), GETV(pi, i),
                   throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pzib(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), size.length(),
                prob.length(), pi.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    size.length(),
    prob.length(),
    pi.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_zib(GETV(x, i), GETV(size, i),
                   GETV(prob, i), GETV(pi, i),
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
NumericVector cpp_qzib(
    const NumericVector& p,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), size.length(),
                prob.length(), pi.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    p.length(),
    size.length(),
    prob.length(),
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
    x[i] = invcdf_zib(GETV(pp, i), GETV(size, i),
                      GETV(prob, i), GETV(pi, i),
                      throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rzib(
    const int& n,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi
  ) {
  
  if (std::min({size.length(), prob.length(), pi.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_zib(GETV(size, i), GETV(prob, i),
                   GETV(pi, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

