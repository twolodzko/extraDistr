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

double pdf_zip(double x, double lambda, double pi) {
  if (lambda <= 0.0 || pi < 0.0 || pi > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(pi))
    return NA_REAL;
  if (x < 0.0 || !isInteger(x) || !R_FINITE(x))
    return 0.0;
  if (x == 0.0)
    return pi + (1.0-pi) * exp(-lambda);
  else
    return (1.0-pi) * R::dpois(x, lambda, false);
}

double cdf_zip(double x, double lambda, double pi) {
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(pi))
    return NA_REAL;
  if (lambda <= 0.0 || pi < 0.0 || pi > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  if (!R_FINITE(x))
    return 1.0;
  return pi + (1.0-pi) * R::ppois(x, lambda, true, false);
}

double invcdf_zip(double p, double lambda, double pi) {
  if (ISNAN(p) || ISNAN(lambda) || ISNAN(pi))
    return NA_REAL;
  if (lambda <= 0.0 || pi < 0.0 || pi > 1.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (p < pi)
    return 0.0;
  else
    return R::qpois((p - pi) / (1.0-pi), lambda, true, false);
}

double rng_zip(double lambda, double pi) {
  if (ISNAN(lambda) || ISNAN(pi) ||
      lambda <= 0.0 || pi < 0.0 || pi > 1.0) {
    Rcpp::warning("NAs produced");
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
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(lambda.length());
  dims.push_back(pi.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_zip(x[i % dims[0]], lambda[i % dims[1]], pi[i % dims[2]]);
  
  if (log_prob)
    p = Rcpp::log(p);
  
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
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(lambda.length());
  dims.push_back(pi.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_zip(x[i % dims[0]], lambda[i % dims[1]], pi[i % dims[2]]);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
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
  
  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(lambda.length());
  dims.push_back(pi.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_zip(pp[i % dims[0]], lambda[i % dims[1]], pi[i % dims[2]]);
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rzip(
    const int& n,
    const NumericVector& lambda,
    const NumericVector& pi
  ) {
  
  std::vector<int> dims;
  dims.push_back(lambda.length());
  dims.push_back(pi.length());
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_zip(lambda[i % dims[0]], pi[i % dims[1]]);
  
  return x;
}

