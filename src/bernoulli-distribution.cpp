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
*  Bernoulli distribution
*
*  Values:
*  x
*
*  Parameters:
*  0 <= p <= 1
*
*/

double pdf_bernoulli(double x, double prob) {
  if (ISNAN(x) || ISNAN(prob))
    return NA_REAL;
  if (prob < 0.0 || prob > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x == 1.0)
    return prob;
  if (x == 0.0)
    return 1.0 - prob;
  Rcpp::warning("improper value of x");
  return 0.0;
}

double cdf_bernoulli(double x, double prob) {
  if (ISNAN(x) || ISNAN(prob))
    return NA_REAL;
  if (prob < 0.0 || prob > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  if (x < 1.0)
    return 1.0 - prob;
  return 1.0;
}

double invcdf_bernoulli(double p, double prob) {
  if (ISNAN(p) || ISNAN(prob))
    return NA_REAL;
  if (prob < 0.0 || prob > 1.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return (p <= 1.0 - prob) ? 0.0 : 1.0;
}

double rng_bernoulli(double p) {
  if (ISNAN(p) || p < 0.0 || p > 1.0) {
    Rcpp::warning("NAs produced");
    return NA_REAL;
  }
  double u = rng_unif();
  return (u > p) ? 0.0 : 1.0;
}


// [[Rcpp::export]]
NumericVector cpp_dbern(
    const NumericVector& x,
    const NumericVector& prob,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(prob.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_bernoulli(x[i % dims[0]], prob[i % dims[1]]);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pbern(
    const NumericVector& x,
    const NumericVector& prob,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(prob.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_bernoulli(x[i % dims[0]], prob[i % dims[1]]);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qbern(
    const NumericVector& p,
    const NumericVector& prob,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(prob.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_bernoulli(pp[i % dims[0]], prob[i % dims[1]]);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rbern(
    const int& n,
    const NumericVector& prob
  ) {
  
  int dims = prob.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_bernoulli(prob[i % dims]);
  
  return x;
}

