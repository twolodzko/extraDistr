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


double pdf_hnorm(double x, double sigma, bool& throw_warning) {
  if (ISNAN(x) || ISNAN(sigma))
    return x+sigma;
  if (sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  return 2.0 * R::dnorm(x, 0.0, sigma, false);
}

double cdf_hnorm(double x, double sigma, bool& throw_warning) {
  if (ISNAN(x) || ISNAN(sigma))
    return x+sigma;
  if (sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  return 2.0 * R::pnorm(x, 0.0, sigma, true, false) - 1.0;
}

double invcdf_hnorm(double p, double sigma, bool& throw_warning) {
  if (ISNAN(p) || ISNAN(sigma))
    return p+sigma;
  if (sigma <= 0.0 || p < 0.0 || p > 1.0) {
    throw_warning = true;
    return NAN;
  }
  return R::qnorm((p+1.0)/2.0, 0.0, sigma, true, false);
}

double rng_hnorm(double sigma, bool& throw_warning) {
  if (ISNAN(sigma) || sigma <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  return abs(R::norm_rand()) * sigma;
}


// [[Rcpp::export]]
NumericVector cpp_dhnorm(
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
    p[i] = pdf_hnorm(x[i % dims[0]], sigma[i % dims[1]],
                     throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_phnorm(
    const NumericVector& x,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(sigma.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_hnorm(x[i % dims[0]], sigma[i % dims[1]],
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
NumericVector cpp_qhnorm(
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
    q[i] = invcdf_hnorm(pp[i % dims[0]], sigma[i % dims[1]],
                        throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rhnorm(
    const int& n,
    const NumericVector& sigma
) {
  
  int dims = sigma.length();
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_hnorm(sigma[i % dims], throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

