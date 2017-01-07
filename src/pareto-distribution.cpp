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
 *  Pareto distribution
 *
 *  Values:
 *  x
 *
 *  Parameters:
 *  a, b > 0
 *
 *  f(x)    = (a*b^a) / x^{a+1}
 *  F(x)    = 1 - (b/x)^a
 *  F^-1(p) = b/(1-p)^{1-a}
 *
 */

double pdf_pareto(double x, double a, double b) {
  if (ISNAN(x) || ISNAN(a) || ISNAN(b))
    return NA_REAL;
  if (a <= 0.0 || b <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < b)
    return 0.0;
  return a * pow(b, a) / pow(x, a+1.0);
}

double logpdf_pareto(double x, double a, double b) {
  if (ISNAN(x) || ISNAN(a) || ISNAN(b))
    return NA_REAL;
  if (a <= 0.0 || b <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < b)
    return R_NegInf;
  return log(a) + log(b)*a - log(x)*(a+1.0);
}

double cdf_pareto(double x, double a, double b) {
  if (ISNAN(x) || ISNAN(a) || ISNAN(b))
    return NA_REAL;
  if (a <= 0.0 || b <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < b)
    return 0.0;
  return 1.0 - pow(b/x, a);
}

double invcdf_pareto(double p, double a, double b) {
  if (ISNAN(p) || ISNAN(a) || ISNAN(b))
    return NA_REAL;
  if (a <= 0.0 || b <= 0.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return b / pow(1.0-p, 1.0/a);
}

double rng_pareto(double a, double b) {
  if (ISNAN(a) || ISNAN(b) || a <= 0.0 || b <= 0.0) {
    Rcpp::warning("NAs produced");
    return NA_REAL;
  }
  double u = rng_unif();
  return b / pow(u, 1.0/a);
}


// [[Rcpp::export]]
NumericVector cpp_dpareto(
    const NumericVector& x,
    const NumericVector& a,
    const NumericVector& b,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_pareto(x[i % dims[0]], a[i % dims[1]], b[i % dims[2]]);

  if (!log_prob)
    p = Rcpp::exp(p);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_ppareto(
    const NumericVector& x,
    const NumericVector& a,
    const NumericVector& b,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_pareto(x[i % dims[0]], a[i % dims[1]], b[i % dims[2]]);

  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qpareto(
    const NumericVector& p,
    const NumericVector& a,
    const NumericVector& b,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(p.length());
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
    x[i] = invcdf_pareto(pp[i % dims[0]], a[i % dims[1]], b[i % dims[2]]);

  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rpareto(
    const int& n,
    const NumericVector& a,
    const NumericVector& b
  ) {

  std::vector<int> dims;
  dims.push_back(a.length());
  dims.push_back(b.length());
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_pareto(a[i % dims[0]], b[i % dims[1]]);

  return x;
}

