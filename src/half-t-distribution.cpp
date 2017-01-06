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
  if (ISNAN(nu) || ISNAN(sigma) || sigma <= 0.0 || nu <= 0.0) {
    Rcpp::warning("NAs produced");
    return NA_REAL;
  }
  return abs(R::rt(nu) * sigma);
}


// [[Rcpp::export]]
NumericVector cpp_dht(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& sigma,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(nu.length());
  dims.push_back(sigma.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++) 
    p[i] = pdf_ht(x[i % dims[0]], nu[i % dims[1]], sigma[i % dims[2]]);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pht(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(nu.length());
  dims.push_back(sigma.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_ht(x[i % dims[0]], nu[i % dims[1]], sigma[i % dims[2]]);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qht(
    const NumericVector& p,
    const NumericVector& nu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(nu.length());
  dims.push_back(sigma.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_ht(pp[i % dims[0]], nu[i % dims[1]], sigma[i % dims[2]]);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rht(
    const int& n,
    const NumericVector& nu,
    const NumericVector& sigma
  ) {
  
  std::vector<int> dims;
  dims.push_back(nu.length());
  dims.push_back(sigma.length());
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_ht(nu[i % dims[0]], sigma[i % dims[1]]);
  
  return x;
}

