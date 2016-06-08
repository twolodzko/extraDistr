#include <Rcpp.h>
#include "namespace.h"
#include "shared.h"


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

double pdf_zib(double x, double n, double p, double pi) {
  if (p < 0 || p > 1 || n < 0 || pi < 0 || pi > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0 || !isInteger(x) || std::isinf(x))
    return 0;
  if (x == 0)
    return pi + (1-pi) * pow(1-p, n);
  else
    return (1-pi) * R::dbinom(x, n, p, false);
}

double cdf_zib(double x, double n, double p, double pi) {
  if (p < 0 || p > 1 || n < 0 || pi < 0 || pi > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0)
    return 0;
  if (std::isinf(x))
    return 1;
  return pi + (1-pi) * R::pbinom(x, n, p, true, false);
}

double invcdf_zib(double pp, double n, double p, double pi) {
  if (p < 0 || p > 1 || n < 0 || pi < 0 || pi > 1 || pp < 0 || pp > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (pp < pi)
    return 0;
  else
    return R::qbinom((pp - pi) / (1-pi), n, p, true, false);
}

double rng_zib(double n, double p, double pi) {
  if (p < 0 || p > 1 || n < 0 || pi < 0 || pi > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double u = R::runif(0, 1);
  if (u < pi)
    return 0;
  else
    return R::rbinom(n, p);
}


// [[Rcpp::export]]
NumericVector cpp_dzib(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi,
    bool log_prob = false
  ) {
  
  int n  = x.length();
  int npi = pi.length();
  int ns = size.length();
  int np = prob.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, npi, ns, np));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_zib(x[i % n], size[i % ns], prob[i % np], pi[i % npi]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pzib(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = x.length();
  int npi = pi.length();
  int ns = size.length();
  int np = prob.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, npi, ns, np));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_zib(x[i % n], size[i % ns], prob[i % np], pi[i % np]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qzib(
    const NumericVector& p,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = p.length();
  int npi = pi.length();
  int ns = size.length();
  int np = prob.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, npi, ns, np));
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1-pp[i];
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_zib(pp[i % n], size[i % ns], prob[i % np], pi[i % np]);
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rzib(
    const int n,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi
  ) {
  
  int npi = pi.length();
  int ns = size.length();
  int np = prob.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_zib(size[i % ns], prob[i % np], pi[i % npi]);
  
  return x;
}

