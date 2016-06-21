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


double pdf_hnorm(double x, double sigma) {
  if (sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  return 2.0 * R::dnorm(x, 0.0, sigma, false);
}

double cdf_hnorm(double x, double sigma) {
  if (sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  return 2.0 * R::pnorm(x, 0.0, sigma, true, false) - 1.0;
}

double invcdf_hnorm(double p, double sigma) {
  if (sigma <= 0.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return R::qnorm((p+1.0)/2.0, 0.0, sigma, true, false);
}

double rng_hnorm(double sigma) {
  if (sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return abs(R::rnorm(0.0, sigma));
}


// [[Rcpp::export]]
NumericVector cpp_dhnorm(
    const NumericVector& x,
    const NumericVector& sigma,
    bool log_prob = false
) {
  
  int n  = x.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, ns));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
      p[i] = pdf_hnorm(x[i % n], sigma[i % ns]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_phnorm(
    const NumericVector& x,
    const NumericVector& sigma,
    bool lower_tail = true, bool log_prob = false
) {
  
  int n  = x.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, ns));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_hnorm(x[i % n], sigma[i % ns]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qhnorm(
    const NumericVector& p,
    const NumericVector& sigma,
    bool lower_tail = true, bool log_prob = false
) {
  
  int n  = p.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, ns));
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1.0 - pp[i];
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_hnorm(pp[i % n], sigma[i % ns]);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rhnorm(
    const int n,
    const NumericVector& sigma
) {
  
  int ns = sigma.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_hnorm(sigma[i % ns]);
  
  return x;
}

