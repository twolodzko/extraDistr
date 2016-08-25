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
 * 
 * x >= 0
 * 
 * Parameters:
 * nu > 0
 * sigma > 0
 * 
 * with nu = 1   returns half-Cauchy
 * with nu = Inf returns half-normal
 * 
 */

double pdf_ht(double x, double nu, double sigma) {
  if (ISNAN(x) || ISNAN(nu) || ISNAN(sigma))
    return NA_REAL;
  if (sigma <= 0.0 || nu <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  return 2.0 * R::dt(x/sigma, nu, false)/sigma;
}

double cdf_ht(double x, double nu, double sigma) {
  if (ISNAN(x) || ISNAN(nu) || ISNAN(sigma))
    return NA_REAL;
  if (sigma <= 0.0 || nu <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  return 2.0 * R::pt(x/sigma, nu, true, false) - 1.0;
}

double invcdf_ht(double p, double nu, double sigma) {
  if (ISNAN(p) || ISNAN(nu) || ISNAN(sigma))
    return NA_REAL;
  if (sigma <= 0.0 || nu <= 0.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return R::qt((p+1.0)/2.0, nu, true, false) * sigma;
}

double rng_ht(double nu, double sigma) {
  if (ISNAN(nu) || ISNAN(sigma))
    return NA_REAL;
  if (sigma <= 0.0 || nu <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return abs(R::rt(nu) * sigma);
}


// [[Rcpp::export]]
NumericVector cpp_dht(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& sigma,
    bool log_prob = false
  ) {
  
  int n  = x.length();
  int nn = nu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, ns));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++) 
    p[i] = pdf_ht(x[i % n], nu[i % nn], sigma[i % ns]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pht(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& sigma,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = x.length();
  int nn = nu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, ns));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_ht(x[i % n], nu[i % nn], sigma[i % ns]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qht(
    const NumericVector& p,
    const NumericVector& nu,
    const NumericVector& sigma,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = p.length();
  int nn = nu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, ns));
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1.0 - pp[i];
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_ht(pp[i % n], nu[i % nn], sigma[i % ns]);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rht(
    const int n,
    const NumericVector& nu,
    const NumericVector& sigma
  ) {
  
  int nn = nu.length();
  int ns = sigma.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_ht(nu[i % nn], sigma[i % ns]);
  
  return x;
}

