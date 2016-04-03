#include <Rcpp.h>
#include "shared.h"
using namespace Rcpp;


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
  if (prob < 0 || prob > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x == 1)
    return prob;
  if (x == 0)
    return 1-prob;
  Rcpp::warning("improper value of x");
  return 0;
}

double cdf_bernoulli(double x, double prob) {
  if (prob < 0 || prob > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0)
    return 0;
  if (x < 1)
    return 1-prob;
  return 1;
}

int invcdf_bernoulli(double p, double prob) {
  if (prob < 0 || prob > 1 || p < 0 || p > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (p <= 1-prob)
    return 0;
  else
    return 1;
}


// [[Rcpp::export]]
NumericVector cpp_dbern(NumericVector x, NumericVector prob,
                        bool log_prob = false) {
  
  int n  = x.length();
  int np = prob.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, np));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_bernoulli(x[i % n], prob[i % np]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pbern(NumericVector x, NumericVector prob,
                        bool lower_tail = true, bool log_prob = false) {
  
  int n  = x.length();
  int np = prob.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, np));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_bernoulli(x[i % n], prob[i % np]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
IntegerVector cpp_qbern(NumericVector p, NumericVector prob,
                        bool lower_tail = true, bool log_prob = false) {
  
  int n  = p.length();
  int np = prob.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, np));
  IntegerVector q(Nmax);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = exp(p[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      p[i] = 1-p[i];
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_bernoulli(p[i % n], prob[i % np]);
  
  return q;
}


// [[Rcpp::export]]
IntegerVector cpp_rbern(int n, NumericVector prob) {
  
  int np = prob.length();
  IntegerVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_bernoulli(prob[i % np]);
  
  return x;
}

