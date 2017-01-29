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
 *  Generalized extreme value distribution
 *
 *  Values:
 *  x
 *
 *  Parameters:
 *  mu
 *  sigma > 0
 *  xi
 *
 *  z = (x-mu)/sigma
 *  where 1+xi*z > 0
 *
 *  f(x)    = { 1/sigma * (1-xi*z)^{-1-1/xi} * exp(-(1-xi*z)^{-1/xi})     if xi != 0
 *            { 1/sigma * exp(-z) * exp(-exp(-z))                         otherwise
 *  F(x)    = { exp(-(1+xi*z)^{1/xi})                                     if xi != 0
 *            { exp(-exp(-z))                                             otherwise
 *  F^-1(p) = { mu - sigma/xi * (1 - (-log(1-p))^xi)                      if xi != 0
 *            { mu - sigma * log(-log(1-p))                               otherwise
 *
 */

inline double pdf_gev(double x, double mu, double sigma,
                      double xi, bool& throw_warning) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(xi))
    return x+mu+sigma+xi;
  if (sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double z = (x-mu)/sigma;
  if (1.0+xi*z > 0.0) {
    if (xi != 0.0)
      return 1.0/sigma * pow(1.0+xi*z, -1.0-(1.0/xi)) * exp(-pow(1.0+xi*z, -1.0/xi));
    else
      return 1.0/sigma * exp(-z) * exp(-exp(-z));
  } else {
    return 0.0;
  }
}

inline double cdf_gev(double x, double mu, double sigma,
                      double xi, bool& throw_warning) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(xi))
    return x+mu+sigma+xi;
  if (sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double z = (x-mu)/sigma;
  if (1.0+xi*z > 0.0) {
    if (xi != 0.0)
      return exp(-pow(1.0+xi*z, -1.0/xi));
    else
      return exp(-exp(-z));
  } else {
    return 0.0;
  }
}

inline double invcdf_gev(double p, double mu, double sigma,
                         double xi, bool& throw_warning) {
  if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma) || ISNAN(xi))
    return p+mu+sigma+xi;
  if (sigma <= 0.0 || !VALID_PROB(p)) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (p == 1.0)
    return R_PosInf;
  if (xi != 0.0)
    return mu - sigma/xi * (1.0 - pow(-log(p), -xi));
  else
    return mu - sigma * log(-log(p));
}

inline double rng_gev(double mu, double sigma, double xi,
                      bool& throw_warning) {
  if (ISNAN(mu) || ISNAN(sigma) || ISNAN(xi) || sigma <= 0.0) {
    Rcpp::warning("NAs produced");
    return NA_REAL;
  }
  double u = rng_unif();
  if (xi != 0.0)
    return mu - sigma/xi * (1.0 - pow(-log(u), -xi));
  else
    return mu - sigma * log(-log(u));
}


// [[Rcpp::export]]
NumericVector cpp_dgev(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& xi,
    const bool& log_prob = false
  ) {

  int Nmax = std::max({
    x.length(),
    mu.length(),
    sigma.length(),
    xi.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_gev(GETV(x, i), GETV(mu, i),
                   GETV(sigma, i), GETV(xi, i),
                   throw_warning);

  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pgev(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& xi,
    bool lower_tail = true, bool log_prob = false
  ) {

  int Nmax = std::max({
    x.length(),
    mu.length(),
    sigma.length(),
    xi.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_gev(GETV(x, i), GETV(mu, i),
                   GETV(sigma, i), GETV(xi, i),
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
NumericVector cpp_qgev(
    const NumericVector& p,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& xi,
    bool lower_tail = true, bool log_prob = false
  ) {

  int Nmax = std::max({
    p.length(),
    mu.length(),
    sigma.length(),
    xi.length()
  });
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;

  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_gev(GETV(pp, i), GETV(mu, i),
                      GETV(sigma, i), GETV(xi, i),
                      throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rgev(
    const int& n,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& xi
  ) {

  NumericVector x(n);

  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_gev(GETV(mu, i), GETV(sigma, i),
                   GETV(xi, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

