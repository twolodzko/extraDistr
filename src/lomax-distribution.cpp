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
*  Lomax distribution
*
*  Values:
*  x > 0
*
*  Parameters:
*  lambda > 0
*  kappa > 0
*
*  f(x)    = lambda*kappa / (1+lambda*x)^(kappa+1)
*  F(x)    = 1-(1+lambda*x)^-kappa
*  F^-1(p) = ((1-p)^(-1/kappa)-1) / lambda
*
*/

double pdf_lomax(double x, double lambda, double kappa) {
  if (lambda <= 0.0 || kappa <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x <= 0.0)
    return 0.0;
  return lambda*kappa / pow(1.0+lambda*x, kappa+1.0);
}

double logpdf_lomax(double x, double lambda, double kappa) {
  if (lambda <= 0.0 || kappa <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x <= 0.0)
    return -INFINITY;
  return log(lambda) + log(kappa) - log(1.0+lambda*x)*(kappa+1.0);
}

double cdf_lomax(double x, double lambda, double kappa) {
  if (lambda <= 0.0 || kappa <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x <= 0.0)
    return 0.0;
  return 1.0 - pow(1.0+lambda*x, -kappa);
}

double invcdf_lomax(double p, double lambda, double kappa) {
  if (lambda <= 0.0 || kappa <= 0.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return (pow(1.0-p, -1.0/kappa)-1.0) / lambda;
}


// [[Rcpp::export]]
NumericVector cpp_dlomax(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& kappa,
    bool log_prob = false
  ) {

  int n = x.length();
  int nl = lambda.length();
  int nk = kappa.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nl, nk));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_lomax(x[i % n], lambda[i % nl], kappa[i % nk]);

  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = exp(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_plomax(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& kappa,
    bool lower_tail = true, bool log_prob = false
  ) {

  int n  = x.length();
  int nl = lambda.length();
  int nk = kappa.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nl, nk));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_lomax(x[i % n], lambda[i % nl], kappa[i % nk]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qlomax(
    const NumericVector& p,
    const NumericVector& lambda,
    const NumericVector& kappa,
    bool lower_tail = true, bool log_prob = false
  ) {

  int n  = p.length();
  int nl = lambda.length();
  int nk = kappa.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nl, nk));
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);

  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1.0 - pp[i];

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_lomax(pp[i % n], lambda[i % nl], kappa[i % nk]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rlomax(
    const int n,
    const NumericVector& lambda,
    const NumericVector& kappa
  ) {

  double u;
  int nl = lambda.length();
  int nk = kappa.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    u = rng_unif();
    x[i] = invcdf_lomax(u, lambda[i % nl], kappa[i % nk]);
  }

  return x;
}

