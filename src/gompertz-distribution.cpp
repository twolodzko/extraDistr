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
*  Gompertz distribution
*
*  Values:
*  x >= 0
*
*  Parameters:
*  a > 0
*  b > 0
*
*  f(x)    = a*exp(b*x - a/b * (exp(bx)-1))
*  F(x)    = 1-exp(-a/b * (exp(bx)-1))
*  F^-1(p) = 1/b * log(1 - b/a * log(1-p))
*
* References:
*
* Lenart, A. (2012). The Gompertz distribution and Maximum Likelihood Estimation
* of its parameters - a revision. MPIDR WORKING PAPER WP 2012-008.
*
*/


double pdf_gompertz(double x, double a, double b) {
  if (ISNAN(x) || ISNAN(a) || ISNAN(b))
    return NAN;
  // if (a <= 0.0 || b <= 0.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }
  if (x < 0.0 || std::isinf(x))
    return 0.0;
  return a * exp(b*x - a/b * (exp(b*x) - 1.0));
}

double cdf_gompertz(double x, double a, double b) {
  if (ISNAN(x) || ISNAN(a) || ISNAN(b))
    return NAN;
  // if (a <= 0.0 || b <= 0.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }
  if (std::isinf(x))
    return 1.0;
  if (x < 0.0)
    return 0.0;
  return 1.0 - exp(-a/b * (exp(b*x) - 1.0));
}

double invcdf_gompertz(double p, double a, double b) {
  if (ISNAN(p) || ISNAN(a) || ISNAN(b))
    return NAN;
  // if (a <= 0.0 || b <= 0.0 || p < 0.0 || p > 1.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }
  return 1.0/b * log(1.0 - b/a * log(1.0-p));
}

double logpdf_gompertz(double x, double a, double b) {
  if (ISNAN(x) || ISNAN(a) || ISNAN(b))
    return NAN;
  // if (a <= 0.0 || b <= 0.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }
  if (x < 0.0 || std::isinf(x))
    return -INFINITY;
  return log(a) + (b*x - a/b * (exp(b*x) - 1.0));
}


// [[Rcpp::export]]
NumericVector cpp_dgompertz(
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
  NumericVector a_n = positive_or_nan(a);
  NumericVector b_n = positive_or_nan(b);

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_gompertz(x[i % n], a_n[i % na], b_n[i % nb]);

  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = exp(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pgompertz(
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
  NumericVector a_n = positive_or_nan(a);
  NumericVector b_n = positive_or_nan(b);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_gompertz(x[i % n], a_n[i % na], b_n[i % nb]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qgompertz(
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
  NumericVector a_n = positive_or_nan(a);
  NumericVector b_n = positive_or_nan(b);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1.0 - pp[i];

  pp = zeroone_or_nan(pp);
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_gompertz(pp[i % n], a_n[i % na], b_n[i % nb]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rgompertz(
    const int n,
    const NumericVector& a,
    const NumericVector& b
  ) {

  double u;
  int na = a.length();
  int nb = b.length();
  NumericVector x(n);
  NumericVector a_n = positive_or_nan(a);
  NumericVector b_n = positive_or_nan(b);

  for (int i = 0; i < n; i++) {
    u = rng_unif();
    x[i] = invcdf_gompertz(u, a_n[i % na], b_n[i % nb]);
  }

  return x;
}

