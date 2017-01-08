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
using std::tan;
using std::atan;


double pdf_hcauchy(double x, double sigma, bool& throw_warning) {
  if (ISNAN(x) || ISNAN(sigma))
    return x+sigma;
  if (sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  return 2.0/(M_PI*(1.0 + pow(x/sigma, 2.0)))/sigma;
}

double cdf_hcauchy(double x, double sigma, bool& throw_warning) {
  if (ISNAN(x) || ISNAN(sigma))
    return x+sigma;
  if (sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  return 2.0/M_PI * atan(x/sigma);
}

double invcdf_hcauchy(double p, double sigma, bool& throw_warning) {
  if (ISNAN(p) || ISNAN(sigma))
    return p+sigma;
  if (sigma <= 0.0 || p < 0.0 || p > 1.0) {
    throw_warning = true;
    return NAN;
  }
  return sigma * tan((M_PI*p)/2.0);
}

double rng_hcauchy(double sigma, bool& throw_warning) {
  if (ISNAN(sigma) || sigma <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  return abs(R::rcauchy(0.0, sigma));
}


// [[Rcpp::export]]
NumericVector cpp_dhcauchy(
    const NumericVector& x,
    const NumericVector& sigma,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(sigma.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_hcauchy(x[i % dims[0]], sigma[i % dims[1]],
                       throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_phcauchy(
    const NumericVector& x,
    const NumericVector& sigma,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(sigma.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_hcauchy(x[i % dims[0]], sigma[i % dims[1]],
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
NumericVector cpp_qhcauchy(
    const NumericVector& p,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(sigma.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_hcauchy(pp[i % dims[0]], sigma[i % dims[1]],
                          throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rhcauchy(
    const int& n,
    const NumericVector& sigma
  ) {
  
  int dims = sigma.length();
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_hcauchy(sigma[i % dims], throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

