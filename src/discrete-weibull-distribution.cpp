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
*  Discrete Weibull distribution
*
*  Values:
*  x >= 0
*
*  Parameters:
*  0 < q < 1
*  beta
*
*  f(x)    = q^x^beta - q^(x+1)^beta
*  F(x)    = 1-q^(x+1)^beta
*  F^-1(p) = ceiling(pow(log(1-p)/log(q), 1/beta) - 1)
*
*  Nakagawa and Osaki (1975), "The Discrete Weibull Distribution",
*  IEEE Transactions on Reliability, R-24, pp. 300-301.
*
*/

double pdf_dweibull(double x, double q, double beta) {
  if (ISNAN(x) || ISNAN(q) || ISNAN(beta))
    return NAN;
  if (q <= 0.0 || q >= 1.0 || beta <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(x) || x < 0.0)
    return 0.0;
  return pow(q, pow(x, beta)) - pow(q, pow(x+1.0, beta));
}

double cdf_dweibull(double x, double q, double beta) {
  if (ISNAN(x) || ISNAN(q) || ISNAN(beta))
    return NAN;
  if (q <= 0.0 || q >= 1.0 || beta <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  return 1.0 - pow(q, pow(x+1.0, beta));
}

double invcdf_dweibull(double p, double q, double beta) {
  if (ISNAN(p) || ISNAN(q) || ISNAN(beta))
    return NA_REAL;
  if (q <= 0.0 || q >= 1.0 || beta <= 0.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (p == 0.0)
    return 0.0;
  return ceil(pow(log(1.0 - p)/log(q), 1.0/beta) - 1.0);
}


// [[Rcpp::export]]
NumericVector cpp_ddweibull(
    const NumericVector& x,
    const NumericVector& q,
    const NumericVector& beta,
    bool log_prob = false
  ) {

  int n  = x.length();
  int nq = q.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nq, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_dweibull(x[i % n], q[i % nq], beta[i % nb]);

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pdweibull(
    const NumericVector& x,
    const NumericVector& q,
    const NumericVector& beta,
    bool lower_tail = true, bool log_prob = false
  ) {

  int n  = x.length();
  int nq = q.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nq, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_dweibull(x[i % n], q[i % nq], beta[i % nb]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qdweibull(
    const NumericVector& p,
    const NumericVector& q,
    const NumericVector& beta,
    bool lower_tail = true, bool log_prob = false
  ) {

  int n  = p.length();
  int nq = q.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nq, nb));
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);

  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1.0 - pp[i];

  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_dweibull(pp[i % n], q[i % nq], beta[i % nb]);

  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rdweibull(
    const int n,
    const NumericVector& q,
    const NumericVector& beta
  ) {

  double u;
  int nq = q.length();
  int nb = beta.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    u = rng_unif();
    x[i] = invcdf_dweibull(u, q[i % nq], beta[i % nb]);
  }

  return x;
}

