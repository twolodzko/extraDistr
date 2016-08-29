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

double pdf_nst(double x, double nu, double mu, double sigma) {
  if (ISNAN(x) || ISNAN(nu) || ISNAN(mu) || ISNAN(sigma))
    return NAN;
  // if (nu <= 0.0 || sigma <= 0.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }
  double z = (x - mu)/sigma;
  return R::dt(z, nu, false)/sigma;
}

double cdf_nst(double x, double nu, double mu, double sigma) {
  if (ISNAN(x) || ISNAN(nu) || ISNAN(mu) || ISNAN(sigma))
    return NAN;
  // if (nu <= 0.0 || sigma <= 0.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }
  double z = (x - mu)/sigma;
  return R::pt(z, nu, true, false);
}

double invcdf_nst(double p, double nu, double mu, double sigma) {
  if (ISNAN(p) || ISNAN(nu) || ISNAN(mu) || ISNAN(sigma))
    return NAN;
  // if (nu <= 0.0 || sigma <= 0.0 || p < 0.0 || p > 1.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }
  return R::qt(p, nu, true, false)*sigma + mu;
}

double rng_nst(double nu, double mu, double sigma) {
  if (ISNAN(nu) || ISNAN(mu) || ISNAN(sigma))
    return NAN;
  // if (nu <= 0.0 || sigma <= 0.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }
  return R::rt(nu)*sigma + mu;
}


// [[Rcpp::export]]
NumericVector cpp_dnst(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma,
    bool log_prob = false
  ) {
  
  int n  = x.length();
  int nn = nu.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, nm, ns));
  NumericVector p(Nmax);
  NumericVector sigma_n = positive_or_nan(sigma);
  NumericVector nu_n = positive_or_nan(nu);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_nst(x[i % n], nu_n[i % nn], mu[i % nm], sigma_n[i % ns]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pnst(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = x.length();
  int nn = nu.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, nm, ns));
  NumericVector p(Nmax);
  NumericVector sigma_n = positive_or_nan(sigma);
  NumericVector nu_n = positive_or_nan(nu);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_nst(x[i % n], nu_n[i % nn], mu[i % nm], sigma_n[i % ns]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qnst(
    const NumericVector& p,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = p.length();
  int nn = nu.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, nm, ns));
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  NumericVector sigma_n = positive_or_nan(sigma);
  NumericVector nu_n = positive_or_nan(nu);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1.0 - pp[i];
  
  pp = zeroone_or_nan(pp);
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_nst(pp[i % n], nu_n[i % nn], mu[i % nm], sigma_n[i % ns]);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rnst(
    const int n,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma
  ) {
  
  int nn = nu.length();
  int nm = mu.length();
  int ns = sigma.length();
  NumericVector x(n);
  NumericVector sigma_n = positive_or_nan(sigma);
  NumericVector nu_n = positive_or_nan(nu);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_nst(nu_n[i % nn], mu[i % nm], sigma_n[i % ns]);
  
  return x;
}

