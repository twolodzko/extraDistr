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
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
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
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
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
    for (int i = 0; i < dims[0]; i++)
      pp[i] = exp(pp[i]);
  
  if (!lower_tail)
    for (int i = 0; i < dims[0]; i++)
      pp[i] = 1.0 - pp[i];
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_bernoulli(pp[i % dims[0]], prob[i % dims[1]]);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rbern(
    const int& n,
    const NumericVector& prob
  ) {
  
  int np = prob.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_bern(prob[i % np]);
  
  return x;
}

