#include <Rcpp.h>
#include "shared.h"
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
*  Non-standard t-distribution
*
*  Values:
*  x
*
*  Parameters:
*  nu > 0
*  mu
*  sigma > 0
*
*/

inline double pdf_nst(double x, double nu, double mu, double sigma,
                      bool& throw_warning) {
  if (ISNAN(x) || ISNAN(nu) || ISNAN(mu) || ISNAN(sigma))
    return x+nu+mu+sigma;
  if (nu <= 0.0 || sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  double z = (x - mu)/sigma;
  return R::dt(z, nu, false)/sigma;
}

inline double cdf_nst(double x, double nu, double mu, double sigma,
                      bool& throw_warning) {
  if (ISNAN(x) || ISNAN(nu) || ISNAN(mu) || ISNAN(sigma))
    return x+nu+mu+sigma;
  if (nu <= 0.0 || sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  double z = (x - mu)/sigma;
  return R::pt(z, nu, true, false);
}

inline double invcdf_nst(double p, double nu, double mu, double sigma,
                         bool& throw_warning) {
  if (ISNAN(p) || ISNAN(nu) || ISNAN(mu) || ISNAN(sigma))
    return p+nu+mu+sigma;
  if (nu <= 0.0 || sigma <= 0.0 || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  return R::qt(p, nu, true, false)*sigma + mu;
}

inline double rng_nst(double nu, double mu, double sigma,
                      bool& throw_warning) {
  if (ISNAN(nu) || ISNAN(mu) || ISNAN(sigma) ||
      nu <= 0.0 || sigma <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  return R::rt(nu)*sigma + mu;
}


// [[Rcpp::export]]
NumericVector cpp_dnst(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& log_prob = false
  ) {
  
  int Nmax = std::max({
    x.length(),
    nu.length(),
    mu.length(),
    sigma.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_nst(GETV(x, i), GETV(nu, i),
                   GETV(mu, i), GETV(sigma, i),
                   throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pnst(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  int Nmax = std::max({
    x.length(),
    nu.length(),
    mu.length(),
    sigma.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_nst(GETV(x, i), GETV(nu, i),
                   GETV(mu, i), GETV(sigma, i),
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
NumericVector cpp_qnst(
    const NumericVector& p,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  int Nmax = std::max({
    p.length(),
    nu.length(),
    mu.length(),
    sigma.length()
  });
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_nst(GETV(pp, i), GETV(nu, i),
                      GETV(mu, i), GETV(sigma, i),
                      throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rnst(
    const int& n,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma
  ) {
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_nst(GETV(nu, i), GETV(mu, i),
                   GETV(sigma, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

