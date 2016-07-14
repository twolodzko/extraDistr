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
*  Kumaraswamy distribution
*
*  Values:
*  x in [0, 1]
*
*  Parameters:
*  a > 0
*  b > 0
*
*  f(x)    = a*b*x^{a-1}*(1-x^a)^{b-1}
*  F(x)    = 1-(1-x^a)^b
*  F^-1(p) = 1-(1-p^{1/b})^{1/a}
*
*/

double pdf_kumar(double x, double a, double b) {
  if (a <= 0.0 || b <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x >= 0.0 && x <= 1.0)
    return a*b * pow(x, a-1.0) * pow(1.0-pow(x, a), b-1.0);
  else
    return 0.0;
}

double cdf_kumar(double x, double a, double b) {
  if (a <= 0.0 || b <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x >= 0.0 && x <= 1.0)
    return 1.0 - pow(1.0 - pow(x, a), b);
  else
    return 0.0;
}

double invcdf_kumar(double p, double a, double b) {
  if (a <= 0.0 || b <= 0.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return pow(1.0 - pow(1.0 - p, 1.0/b), 1.0/a);
}

double logpdf_kumar(double x, double a, double b) {
  if (a <= 0.0 || b <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x >= 0.0 && x <= 1.0)
    return log(a) + log(b) + log(x)*(a-1.0) + log(1.0 - pow(x, a))*(b-1.0);
  else
    return -INFINITY;
}


// [[Rcpp::export]]
NumericVector cpp_dkumar(
    const NumericVector& x,
    const NumericVector& a,
    const NumericVector& b,
    bool log_prob = false
  ) {

  int n  = x.length();
  int na = a.length();
  int nb = b.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_kumar(x[i % n], a[i % na], b[i % nb]);

  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = exp(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pkumar(
    const NumericVector& x,
    const NumericVector& a,
    const NumericVector& b,
    bool lower_tail = true, bool log_prob = false
  ) {

  int n  = x.length();
  int na = a.length();
  int nb = b.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_kumar(x[i % n], a[i % na], b[i % nb]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qkumar(
    const NumericVector& p,
    const NumericVector& a,
    const NumericVector& b,
    bool lower_tail = true, bool log_prob = false
  ) {

  int n  = p.length();
  int na = a.length();
  int nb = b.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);

  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1.0 - pp[i];

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_kumar(pp[i % n], a[i % na], b[i % nb]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rkumar(
    const int n,
    const NumericVector& a,
    const NumericVector& b
  ) {

  double u;
  int na = a.length();
  int nb = b.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    u = rng_unif();
    x[i] = invcdf_kumar(u, a[i % na], b[i % nb]);
  }

  return x;
}

