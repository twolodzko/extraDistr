#include <Rcpp.h>
#include "shared.h"
using namespace Rcpp;


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
  if (q <= 0 || q >= 1 || beta <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(x) || x < 0)
    return 0;
  return pow(q, pow(x, beta)) - pow(q, pow(x+1, beta));
}

double cdf_dweibull(double x, double q, double beta) {
  if (q <= 0 || q >= 1 || beta <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0)
    return 0;
  return 1-pow(q, pow(x+1, beta));
}

double invcdf_dweibull(double p, double q, double beta) {
  if (q <= 0 || q >= 1 || beta <= 0 || p < 0 || p > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (p == 0)
    return 0;
  return ceil(pow(log(1-p)/log(q), 1/beta) - 1);
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
      p[i] = 1-p[i];

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
      pp[i] = 1-pp[i];

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
    u = R::runif(0, 1);
    x[i] = invcdf_dweibull(u, q[i % nq], beta[i % nb]);
  }

  return x;
}

