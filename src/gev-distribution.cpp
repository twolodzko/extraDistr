#include <Rcpp.h>
#include "shared.h"

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using std::sin;
using std::cos;
using std::tan;
using std::atan;
using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;


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

double pdf_gev(double x, double mu, double sigma, double xi) {
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

double cdf_gev(double x, double mu, double sigma, double xi) {
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

double invcdf_gev(double p, double mu, double sigma, double xi) {
  if (sigma <= 0.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (p == 1.0)
    return INFINITY;
  if (xi != 0.0)
    return mu - sigma/xi * (1.0 - pow(-log(p), -xi));
  else
    return mu - sigma * log(-log(p));
}


// [[Rcpp::export]]
NumericVector cpp_dgev(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& xi,
    bool log_prob = false
  ) {

  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int nx = xi.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, nx));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_gev(x[i % n], mu[i % nm], sigma[i % ns], xi[i % nx]);

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

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

  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int nx = xi.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, nx));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_gev(x[i % n], mu[i % nm], sigma[i % ns], xi[i % nx]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

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

  int n  = p.length();
  int nm = mu.length();
  int ns = sigma.length();
  int nx = xi.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, nx));
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);

  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1.0 - pp[i];

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_gev(pp[i % n], mu[i % nm], sigma[i % ns], xi[i % nx]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rgev(
    const int n,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& xi
  ) {

  double u;
  int nm = mu.length();
  int ns = sigma.length();
  int nx = xi.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    u = rng_unif();
    x[i] = invcdf_gev(u, mu[i % nm], sigma[i % ns], xi[i % nx]);
  }

  return x;
}

