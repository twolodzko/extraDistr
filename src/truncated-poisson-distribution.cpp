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
    return NAN;
  // if (lambda < 0.0 || s < 0.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }
  if (!isInteger(x) || x < 0.0 || std::isinf(x))
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
    return NAN;
  // if (lambda <= 0.0 || s < 0.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }
  
  if (x < 0.0)
    return 0.0;
  if (x == INFINITY)
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
    return NAN;
  // if (lambda < 0.0 || s < 0.0 || p < 0.0 || p > 1.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }

  double z;
  
  if (s == 0.0) {
    if (p == 0.0)
      return 0.0;
    if (p == 1.0)
      return INFINITY;
    
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
    return NAN;
  // if (lambda < 0.0 || s < 0.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }

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
    bool log_prob = false
) {
  
  int n  = x.length();
  int nl = lambda.length();
  int ns = s.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nl, ns));
  NumericVector p(Nmax);
  NumericVector lambda_n = nonneg_or_nan(lambda);
  NumericVector s_n = nonneg_or_nan(s);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_tpois(x[i % n], lambda_n[i % nl], s_n[i % ns]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_ptpois(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& s,
    bool lower_tail = true, bool log_prob = false
) {
  
  int n  = x.length();
  int nl = lambda.length();
  int ns = s.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nl, ns));
  NumericVector p(Nmax);
  NumericVector lambda_n = nonneg_or_nan(lambda);
  NumericVector s_n = nonneg_or_nan(s);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_tpois(x[i % n], lambda_n[i % nl], s_n[i % ns]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qtpois(
    const NumericVector& p,
    const NumericVector& lambda,
    const NumericVector& s,
    bool lower_tail = true, bool log_prob = false
) {
  
  int n  = p.length();
  int nl = lambda.length();
  int ns = s.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nl, ns));
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  NumericVector lambda_n = nonneg_or_nan(lambda);
  NumericVector s_n = nonneg_or_nan(s);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1.0 - pp[i];
  
  pp = zeroone_or_nan(pp);
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_tpois(pp[i % n], lambda_n[i % nl], s_n[i % ns]);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rtpois(
    const int n,
    const NumericVector& lambda,
    const NumericVector& s
) {
  
  int nl = lambda.length();
  int ns = s.length();
  NumericVector x(n);
  NumericVector lambda_n = nonneg_or_nan(lambda);
  NumericVector s_n = nonneg_or_nan(s);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_tpois(lambda_n[i % nl], s_n[i % ns]);
  
  return x;
}

