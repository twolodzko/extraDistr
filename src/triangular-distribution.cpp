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
*  Triangular distribution
*
*  Values:
*  x
*
*  Parameters:
*  a
*  b > a
*  a <= c <= b
*
*  f(x)    = { (2*(x-a)) / ((b-a)*(c-a))  x < c
*            { 2/(b-a)                    x = c
*            { (2*(b-x)) / ((b-a)*(b-c))  x > c
*  F(x)    = { (x-a)^2 / ((b-a)*(c-a))
*            { 1 - ((b-x)^2 / ((b-a)*(b-c)))
*  F^-1(p) = { a + sqrt(p*(b-a)*(c-a))    p < (c-a)/(b-a)
*            { b - sqrt((1-p)*(b-a)*(b-c));
*/

inline double pdf_triangular(double x, double a, double b,
                             double c, bool& throw_warning) {
  if (ISNAN(x) || ISNAN(a) || ISNAN(b) || ISNAN(c))
    return x+a+b+c;
  if (a > c || c > b || a == b) {
    throw_warning = true;
    return NAN;
  }
  if (x < a || x > b) {
    return 0.0;
  } else if (x < c) {
    return 2.0*(x-a) / ((b-a)*(c-a));
  } else if (x > c) {
    return 2.0*(b-x) / ((b-a)*(b-c));
  } else {
    return 2.0/(b-a);
  }
}

inline double cdf_triangular(double x, double a, double b,
                             double c, bool& throw_warning) {
  if (ISNAN(x) || ISNAN(a) || ISNAN(b) || ISNAN(c))
    return x+a+b+c;
  if (a > c || c > b || a == b) {
    throw_warning = true;
    return NAN;
  }
  if (x < a) {
    return 0.0;
  } else if (x >= b) {
    return 1.0;
  } else if (x <= c) {
    return pow(x-a, 2.0) / ((b-a)*(c-a));
  } else {
    return 1.0 - (pow(b-x, 2.0) / ((b-a)*(b-c)));
  }
}

inline double invcdf_triangular(double p, double a, double b,
                                double c, bool& throw_warning) {
  if (ISNAN(p) || ISNAN(a) || ISNAN(b) || ISNAN(c))
    return p+a+b+c;
  if (a > c || c > b || a == b || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  double fc = (c-a)/(b-a);
  if (p < fc)
    return a + sqrt(p*(b-a)*(c-a));
  return b - sqrt((1.0-p)*(b-a)*(b-c));
}

inline double rng_triangular(double a, double b, double c,
                             bool& throw_warning) {
  if (ISNAN(a) || ISNAN(b) || ISNAN(c) ||
      a > c || c > b || a == b) {
    throw_warning = true;
    return NA_REAL;
  }
  double u, v, r, cc;
  r = b - a;
  cc = (c-a)/r;
  u = rng_unif();
  v = rng_unif();
  return ((1.0-cc) * std::min(u, v) + cc * std::max(u, v)) * r + a;
}


// [[Rcpp::export]]
NumericVector cpp_dtriang(
    const NumericVector& x,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  dims.push_back(c.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_triangular(x[i % dims[0]], a[i % dims[1]],
                          b[i % dims[2]], c[i % dims[3]],
                          throw_warning);

  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_ptriang(
    const NumericVector& x,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  dims.push_back(c.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_triangular(x[i % dims[0]], a[i % dims[1]],
                          b[i % dims[2]], c[i % dims[3]],
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
NumericVector cpp_qtriang(
    const NumericVector& p,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  dims.push_back(c.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;

  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;

  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_triangular(pp[i % dims[0]], a[i % dims[1]],
                             b[i % dims[2]], c[i % dims[3]],
                             throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rtriang(
    const int& n,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c
  ) {

  std::vector<int> dims;
  dims.push_back(a.length());
  dims.push_back(b.length());
  dims.push_back(c.length());
  NumericVector x(n);
  
  bool throw_warning = false;

  for (int i = 0; i < n; i++)
    x[i] = rng_triangular(a[i % dims[0]], b[i % dims[1]],
                          c[i % dims[2]], throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

