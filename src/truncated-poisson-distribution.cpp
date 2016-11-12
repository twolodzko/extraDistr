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


double pdf_tpois(double x, double lambda, double s) {
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(s))
    return NA_REAL;
  if (lambda < 0.0 || s < 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(x) || x < 0.0 || !R_FINITE(x))
    return 0.0;
  
  if (s == 0.0 && x <= s)
    return 0.0;
  if (s > 0.0 && x > s)
    return 0.0;
  if (s == 0.0)
    return pow(lambda, x) / (factorial(x) * (exp(lambda) - 1.0));
  return R::dpois(x, lambda, false) / R::ppois(s, lambda, true, false);
}

double cdf_tpois(double x, double lambda, double s) {
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(s))
    return NA_REAL;
  if (lambda <= 0.0 || s < 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (x < 0.0)
    return 0.0;
  if (!R_FINITE(x))
    return 1.0;
  
  if (s == 0.0 && x <= s)
    return 0.0;
  if (s > 0.0 && x > s)
    return 1.0;
  if (s == 0.0)
    return R::ppois(x, lambda, true, false) / (1.0 - exp(-lambda));
  return R::ppois(x, lambda, true, false) / R::ppois(s, lambda, true, false);
}

double invcdf_tpois(double p, double lambda, double s) {
  if (ISNAN(p) || ISNAN(lambda) || ISNAN(s))
    return NA_REAL;
  if (lambda < 0.0 || s < 0.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }

  double z;
  
  if (s == 0.0) {
    if (p == 0.0)
      return 0.0;
    if (p == 1.0)
      return R_PosInf;
    
    z = exp(-lambda);
    return R::qpois(p*(1.0-z) + z, lambda, true, false);
  }

  if (p == 0.0)
    return 0.0;
  if (p == 1.0)
    return s;
  
  z = R::ppois(s, lambda, true, false);
  return R::qpois(p*z, lambda, true, false);
}

double rng_tpois(double lambda, double s) {
  if (ISNAN(lambda) || ISNAN(s))
    return NA_REAL;
  if (lambda < 0.0 || s < 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }

  double z, u;
  
  if (s == 0.0) {
    z = exp(-lambda);
    u = R::runif(z, 1.0);
    return R::qpois(u, lambda, true, false);
  }
  
  z = R::ppois(s, lambda, true, false);
  u = R::runif(0.0, z);
  return R::qpois(u, lambda, true, false);
}


// [[Rcpp::export]]
NumericVector cpp_dtpois(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& s,
    const bool& log_prob = false
) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(lambda.length());
  dims.push_back(s.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_tpois(x[i % dims[0]], lambda[i % dims[1]], s[i % dims[2]]);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_ptpois(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& s,
    const bool& lower_tail = true,
    const bool& log_prob = false
) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(lambda.length());
  dims.push_back(s.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_tpois(x[i % dims[0]], lambda[i % dims[1]], s[i % dims[2]]);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qtpois(
    const NumericVector& p,
    const NumericVector& lambda,
    const NumericVector& s,
    const bool& lower_tail = true,
    const bool& log_prob = false
) {
  
  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(lambda.length());
  dims.push_back(s.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_tpois(pp[i % dims[0]], lambda[i % dims[1]], s[i % dims[2]]);
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rtpois(
    const int& n,
    const NumericVector& lambda,
    const NumericVector& s
) {
  
  std::vector<int> dims;
  dims.push_back(lambda.length());
  dims.push_back(s.length());
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_tpois(lambda[i % dims[0]], s[i % dims[1]]);
  
  return x;
}

