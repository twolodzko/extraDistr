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


double pdf_tbinom(double x, double size, double prob, double a, double b) {
  if (ISNAN(x) || ISNAN(size) || ISNAN(prob) || ISNAN(a) || ISNAN(b))
    return NA_REAL;
  if (size < 0.0 || prob < 0.0 || prob > 1.0 || b < a || !isInteger(size, false)) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (!isInteger(x) || x < 0.0 || x <= a || x > b || x > size)
    return 0.0;
  
  double pa, pb;
  pa = R::pbinom(a, size, prob, true, false);
  pb = R::pbinom(b, size, prob, true, false);
  
  return R::dbinom(x, size, prob, false) / (pb-pa);
}

double cdf_tbinom(double x, double size, double prob, double a, double b) {
  if (ISNAN(x) || ISNAN(size) || ISNAN(prob) || ISNAN(a) || ISNAN(b))
    return NA_REAL;
  if (size < 0.0 || prob < 0.0 || prob > 1.0 || b < a || !isInteger(size, false)) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (x < 0.0 || x <= a)
    return 0.0;
  if (x > b || x >= size)
    return 1.0;
  
  double pa, pb;
  pa = R::pbinom(a, size, prob, true, false);
  pb = R::pbinom(b, size, prob, true, false);
  
  return (R::pbinom(x, size, prob, true, false) - pa) / (pb-pa);
}

double invcdf_tbinom(double p, double size, double prob, double a, double b) {
  if (ISNAN(p) || ISNAN(size) || ISNAN(prob) || ISNAN(a) || ISNAN(b))
    return NA_REAL;
  if (size < 0.0 || prob < 0.0 || prob > 1.0 || b < a ||
      !isInteger(size, false) || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (p == 0.0)
    return std::max(a, 0.0);
  if (p == 1.0)
    return std::min(size, b);
  
  double pa, pb;
  pa = R::pbinom(a, size, prob, true, false);
  pb = R::pbinom(b, size, prob, true, false);
  
  return R::qbinom(pa + p*(pb-pa), size, prob, true, false);
}

double rng_tbinom(double size, double prob, double a, double b) {
  if (ISNAN(size) || ISNAN(prob) || ISNAN(a) || ISNAN(b) ||
      size < 0.0 || prob < 0.0 || prob > 1.0 || b < a || !isInteger(size, false)) {
    Rcpp::warning("NAs produced");
    return NA_REAL;
  }
  
  double u, pa, pb;
  pa = R::pbinom(a, size, prob, true, false);
  pb = R::pbinom(b, size, prob, true, false);
  
  u = R::runif(pa, pb);
  return R::qbinom(u, size, prob, true, false);
}


// [[Rcpp::export]]
NumericVector cpp_dtbinom(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& a,
    const NumericVector& b,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(size.length());
  dims.push_back(prob.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_tbinom(x[i % dims[0]], size[i % dims[1]], prob[i % dims[2]],
                      a[i % dims[3]], b[i % dims[4]]);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_ptbinom(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& a,
    const NumericVector& b,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(size.length());
  dims.push_back(prob.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_tbinom(x[i % dims[0]], size[i % dims[1]], prob[i % dims[2]],
                      a[i % dims[3]], b[i % dims[4]]);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qtbinom(
    const NumericVector& p,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& a,
    const NumericVector& b,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(size.length());
  dims.push_back(prob.length());
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
    x[i] = invcdf_tbinom(pp[i % dims[0]], size[i % dims[1]], prob[i % dims[2]],
                         a[i % dims[3]], b[i % dims[4]]);
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rtbinom(
    const int& n,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& a,
    const NumericVector& b
  ) {
  
  std::vector<int> dims;
  dims.push_back(size.length());
  dims.push_back(prob.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_tbinom(size[i % dims[0]], prob[i % dims[1]], a[i % dims[2]], b[i % dims[3]]);
  
  return x;
}

