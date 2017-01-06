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


double pdf_tpois(double x, double lambda, double a, double b) {
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(a) || ISNAN(b))
    return NA_REAL;
  if (lambda < 0.0 || b < a) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (!isInteger(x) || x < 0.0 || x <= a || x > b || !R_FINITE(x))
    return 0.0;
  
  // if (a == 0.0 && b == R_PosInf)
  //   return pow(lambda, x) / (factorial(x) * (exp(lambda) - 1.0));
  
  double pa, pb;
  pa = R::ppois(a, lambda, true, false);
  pb = R::ppois(b, lambda, true, false);
  
  return R::dpois(x, lambda, false) / (pb-pa);
}

double cdf_tpois(double x, double lambda, double a, double b) {
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(a) || ISNAN(b))
    return NA_REAL;
  if (lambda <= 0.0 || b < a) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (x < 0.0 || x <= a)
    return 0.0;
  if (x > b || !R_FINITE(x))
    return 1.0;
  
  // if (a == 0.0 && b == R_PosInf)
  //   return R::ppois(x, lambda, true, false) / (1.0 - exp(-lambda));
  
  double pa, pb;
  pa = R::ppois(a, lambda, true, false);
  pb = R::ppois(b, lambda, true, false);

  return (R::ppois(x, lambda, true, false) - pa) / (pb-pa);
}

double invcdf_tpois(double p, double lambda, double a, double b) {
  if (ISNAN(p) || ISNAN(lambda) || ISNAN(a) || ISNAN(b))
    return NA_REAL;
  if (lambda < 0.0 || b < a || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }

  if (p == 0.0)
    return std::max(a, 0.0);
  if (p == 1.0)
    return b;
  
  double pa, pb;
  pa = R::ppois(a, lambda, true, false);
  pb = R::ppois(b, lambda, true, false);
  
  return R::qpois(pa + p*(pb-pa), lambda, true, false);
}

double rng_tpois(double lambda, double a, double b) {
  if (ISNAN(lambda) || ISNAN(a) || ISNAN(b) ||
      lambda < 0.0 || b < a) {
    Rcpp::warning("NAs produced");
    return NA_REAL;
  }

  double u, pa, pb;
  pa = R::ppois(a, lambda, true, false);
  pb = R::ppois(b, lambda, true, false);
  
  u = R::runif(pa, pb);
  return R::qpois(u, lambda, true, false);
}


// [[Rcpp::export]]
NumericVector cpp_dtpois(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& a,
    const NumericVector& b,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(lambda.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_tpois(x[i % dims[0]], lambda[i % dims[1]], a[i % dims[2]], b[i % dims[3]]);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_ptpois(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& a,
    const NumericVector& b,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(lambda.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_tpois(x[i % dims[0]], lambda[i % dims[1]], a[i % dims[2]], b[i % dims[3]]);
  
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
    const NumericVector& a,
    const NumericVector& b,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(lambda.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_tpois(pp[i % dims[0]], lambda[i % dims[1]], a[i % dims[2]], b[i % dims[3]]);
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rtpois(
    const int& n,
    const NumericVector& lambda,
    const NumericVector& a,
    const NumericVector& b
  ) {
  
  std::vector<int> dims;
  dims.push_back(lambda.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_tpois(lambda[i % dims[0]], a[i % dims[1]], b[i % dims[2]]);
  
  return x;
}

