#include <Rcpp.h>
#include "shared.h"
using namespace Rcpp;

/*
* Zero-inflated Poisson distribution
* 
* Parameters:
* lambda > 0
* 0 <= pi <= 1
* 
* Values:
* x >= 0
*
*/

double pdf_zip(double x, double lambda, double pi) {
  if (lambda <= 0 || pi < 0 || pi > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0 || !isInteger(x) || std::isinf(x))
    return 0;
  if (x == 0)
    return pi + (1-pi) * exp(-lambda);
  else
    return (1-pi) * R::dpois(x, lambda, false);
}

double cdf_zip(double x, double lambda, double pi) {
  if (lambda <= 0 || pi < 0 || pi > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0)
    return 0;
  if (std::isinf(x))
    return 1;
  return pi + (1-pi) * R::ppois(x, lambda, true, false);
}

double invcdf_zip(double p, double lambda, double pi) {
  if (lambda <= 0 || pi < 0 || pi > 1 || p < 0 || p > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (p < pi)
    return 0;
  else
    return R::qpois((p - pi) / (1-pi), lambda, true, false);
}

double rng_zip(double lambda, double pi) {
  if (lambda <= 0 || pi < 0 || pi > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double u = R::runif(0, 1);
  if (u < pi)
    return 0;
  else
    return R::rpois(lambda);
}


// [[Rcpp::export]]
NumericVector cpp_dzip(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& pi,
    bool log_prob = false
  ) {
  
  int n  = x.length();
  int np = pi.length();
  int nl = lambda.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, np, nl));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_zip(x[i % n], lambda[i % nl], pi[i % np]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pzip(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& pi,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = x.length();
  int np = pi.length();
  int nl = lambda.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, np, nl));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_zip(x[i % n], lambda[i % nl], pi[i % np]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qzip(
    const NumericVector& p,
    const NumericVector& lambda,
    const NumericVector& pi,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = p.length();
  int np = pi.length();
  int nl = lambda.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, np, nl));
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1-pp[i];
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_zip(pp[i % n], lambda[i % nl], pi[i % np]);
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rzip(
    const int n,
    const NumericVector& lambda,
    const NumericVector& pi
  ) {
  
  int np = pi.length();
  int nl = lambda.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_zip(lambda[i % nl], pi[i % np]);
  
  return x;
}

