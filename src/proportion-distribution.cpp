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
*  Re-parametrized beta distribution
*
*  Values:
*  x
*
*  Parameters:
*  0 <= mean <= 1
*  size > 0
*
*/

inline double pdf_prop(double x, double size, double mean,
                       bool& throw_warning) {
  if (ISNAN(x) || ISNAN(size) || ISNAN(mean))
    return x+size+mean;
  if (size <= 0.0 || mean <= 0.0 || mean >= 1.0) {
    throw_warning = true;
    return NAN;
  }
  return R::dbeta(x, size*mean+1.0, size*(1.0-mean)+1.0, false);
}

inline double cdf_prop(double x, double size, double mean,
                       bool& throw_warning) {
  if (ISNAN(x) || ISNAN(size) || ISNAN(mean))
    return x+size+mean;
  if (size <= 0.0 || mean <= 0.0 || mean >= 1.0) {
    throw_warning = true;
    return NAN;
  }
  return R::pbeta(x, size*mean+1.0, size*(1.0-mean)+1.0, true, false);
}

inline double invcdf_prop(double p, double size, double mean,
                          bool& throw_warning) {
  if (ISNAN(p) || ISNAN(size) || ISNAN(mean))
    return p+size+mean;
  if (size <= 0.0 || mean <= 0.0 || mean >= 1.0 || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  return R::qbeta(p, size*mean+1.0, size*(1.0-mean)+1.0, true, false);
}

inline double rng_prop(double size, double mean,
                       bool& throw_warning) {
  if (ISNAN(size) || ISNAN(mean) ||
      size <= 0.0 || mean <= 0.0 || mean >= 1.0) {
    throw_warning = true;
    return NA_REAL;
  }
  return R::rbeta(size*mean+1.0, size*(1.0-mean)+1.0);
}


// [[Rcpp::export]]
NumericVector cpp_dprop(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& mean,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mean.length());
  dims.push_back(size.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_prop(GETV(x, i), GETV(size, i),
                    GETV(mean, i), throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pprop(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& mean,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mean.length());
  dims.push_back(size.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_prop(GETV(x, i), GETV(size, i),
                    GETV(mean, i), throw_warning);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qprop(
    const NumericVector& p,
    const NumericVector& size,
    const NumericVector& mean,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(mean.length());
  dims.push_back(size.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_prop(GETV(pp, i), GETV(size, i),
                       GETV(mean, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rprop(
    const int& n,
    const NumericVector& size,
    const NumericVector& mean
  ) {
  
  std::vector<int> dims;
  dims.push_back(mean.length());
  dims.push_back(size.length());
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_prop(GETV(size, i), GETV(mean, i),
                    throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

