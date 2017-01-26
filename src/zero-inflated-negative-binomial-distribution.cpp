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

inline double pdf_zinb(double x, double r, double p, double pi,
                       bool& throw_warning) {
  if (ISNAN(x) || ISNAN(r) || ISNAN(p) || ISNAN(pi))
    return x+r+p+pi;
  if (!VALID_PROB(p) || r < 0.0 || !VALID_PROB(pi) || !isInteger(r, false)) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0 || !isInteger(x) || !R_FINITE(x))
    return 0.0;
  if (x == 0.0)
    return pi + (1.0-pi) * pow(p, r);
  else
    return (1.0-pi) * R::dnbinom(x, r, p, false);
}

inline double cdf_zinb(double x, double r, double p, double pi,
                       bool& throw_warning) {
  if (ISNAN(x) || ISNAN(r) || ISNAN(p) || ISNAN(pi))
    return x+r+p+pi;
  if (!VALID_PROB(p) || r < 0.0 || !VALID_PROB(pi) || !isInteger(r, false)) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  if (!R_FINITE(x))
    return 1.0;
  return pi + (1.0-pi) * R::pnbinom(x, r, p, true, false);
}

inline double invcdf_zinb(double pp, double r, double p, double pi,
                          bool& throw_warning) {
  if (ISNAN(pp) || ISNAN(r) || ISNAN(p) || ISNAN(pi))
    return pp+r+p+pi;
  if (!VALID_PROB(p) || r < 0.0 || !VALID_PROB(pi) ||
      !isInteger(r, false) || !VALID_PROB(pp)) {
    throw_warning = true;
    return NAN;
  }
  if (pp < pi)
    return 0.0;
  else
    return R::qnbinom((pp - pi) / (1.0-pi), r, p, true, false);
}

inline double rng_zinb(double r, double p, double pi,
                       bool& throw_warning) {
  if (ISNAN(r) || ISNAN(p) || ISNAN(pi) || !VALID_PROB(p) ||
      r < 0.0 || !VALID_PROB(pi) || !isInteger(r, false)) {
    throw_warning = true;
    return NA_REAL;
  }
  double u = rng_unif();
  if (u < pi)
    return 0.0;
  else
    return R::rnbinom(r, p);
}


// [[Rcpp::export]]
NumericVector cpp_dzinb(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi,
    const bool& log_prob = false
  ) {
  
  int Nmax = std::max({
    x.length(),
    size.length(),
    prob.length(),
    pi.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_zinb(GETV(x, i), GETV(size, i),
                    GETV(prob, i), GETV(pi, i),
                    throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pzinb(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  int Nmax = std::max({
    x.length(),
    size.length(),
    prob.length(),
    pi.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_zinb(GETV(x, i), GETV(size, i),
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
NumericVector cpp_qzinb(
    const NumericVector& p,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
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
    x[i] = invcdf_zinb(GETV(pp, i), GETV(size, i),
                       GETV(prob, i), GETV(pi, i),
                       throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rzinb(
    const int& n,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi
  ) {
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_zinb(GETV(size, i), GETV(prob, i),
                    GETV(pi, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

