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

double pdf_zinb(double x, double r, double p, double pi) {
  if (p < 0 || p > 1 || r < 0 || pi < 0 || pi > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0 || !isInteger(x))
    return 0;
  if (x == 0)
    return pi + (1-pi) * pow(p, r);
  else
    return (1-pi) * R::dnbinom(x, r, p, false);
}

double cdf_zinb(double x, double r, double p, double pi) {
  if (p < 0 || p > 1 || r < 0 || pi < 0 || pi > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0)
    return 0;
  else
    return pi + (1-pi) * R::pnbinom(x, r, p, true, false);
}

double invcdf_zinb(double pp, double r, double p, double pi) {
  if (p < 0 || p > 1 || r < 0 || pi < 0 || pi > 1 || pp < 0 || pp > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (p < pi)
    return 0;
  else
    return R::qnbinom((pp - pi) / (1-pi), r, p, true, false);
}

double rng_zinb(double r, double p, double pi) {
  if (p < 0 || p > 1 || r < 0 || pi < 0 || pi > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double u = R::runif(0, 1);
  if (u < pi)
    return 0;
  else
    return R::rnbinom(r, p);
}


// [[Rcpp::export]]
NumericVector cpp_dzinb(NumericVector x,
                        NumericVector size, NumericVector prob, NumericVector pi,
                        bool log_prob = false) {
  
  int n  = x.length();
  int npi = pi.length();
  int ns = size.length();
  int np = prob.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, npi, ns, np));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_zinb(x[i % n], size[i % ns], prob[i % np], pi[i % npi]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pzinb(NumericVector x,
                        NumericVector size, NumericVector prob, NumericVector pi,
                        bool lower_tail = true, bool log_prob = false) {
  
  int n  = x.length();
  int npi = pi.length();
  int ns = size.length();
  int np = prob.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, npi, ns, np));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_zinb(x[i % n], size[i % ns], prob[i % np], pi[i % np]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qzinb(NumericVector p,
                       NumericVector size, NumericVector prob, NumericVector pi,
                       bool lower_tail = true, bool log_prob = false) {
  
  int n  = p.length();
  int npi = pi.length();
  int ns = size.length();
  int np = prob.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, npi, ns, np));
  NumericVector x(Nmax);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = exp(p[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      p[i] = 1-p[i];
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_zinb(p[i % n], size[i % ns], prob[i % np], pi[i % np]);
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rzinb(int n,
                        NumericVector size, NumericVector prob, NumericVector pi) {
  
  int npi = pi.length();
  int ns = size.length();
  int np = prob.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_zinb(size[i % ns], prob[i % np], pi[i % npi]);
  
  return x;
}

