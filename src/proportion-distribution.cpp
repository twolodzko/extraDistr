#include <Rcpp.h>
#include "namespace.h"


/*
*  Re-parametrized beta distribution
*
*  Values:
*  x
*
*  Parameters:
*  0 <= mean <= 1
*  size > 0
*
*/

double pdf_prop(double x, double size, double mean, bool log_p) {
  if (size <= 0 || mean < 0 || mean > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return R::dbeta(x, size*mean+1, size*(1-mean)+1, log_p);
}

double cdf_prop(double x, double size, double mean, bool lower_tail, bool log_p) {
  if (size <= 0 || mean < 0 || mean > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return R::pbeta(x, size*mean+1, size*(1-mean)+1, lower_tail, log_p);
}

double invcdf_prop(double p, double size, double mean, bool lower_tail, bool log_p) {
  if (size <= 0 || mean < 0 || mean > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return R::qbeta(p, size*mean+1, size*(1-mean)+1, lower_tail, log_p);
}

double rng_prop(double size, double mean) {
  if (size <= 0 || mean < 0 || mean > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return R::rbeta(size*mean+1, size*(1-mean)+1);
}


// [[Rcpp::export]]
NumericVector cpp_dprop(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& mean,
    bool log_prob = false
  ) {
  
  int n  = x.length();
  int nm = mean.length();
  int ns = size.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_prop(x[i % n], size[i % ns], mean[i % nm], log_prob);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pprop(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& mean,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = x.length();
  int nm = mean.length();
  int ns = size.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_prop(x[i % n], size[i % ns], mean[i % nm], lower_tail, log_prob);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qprop(
    const NumericVector& p,
    const NumericVector& size,
    const NumericVector& mean,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = p.length();
  int nm = mean.length();
  int ns = size.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns));
  NumericVector q(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_prop(p[i % n], size[i % ns], mean[i % nm], lower_tail, log_prob);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rprop(
    const int n,
    const NumericVector& size,
    const NumericVector& mean
  ) {
  
  int nm = mean.length();
  int ns = size.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_prop(size[i % ns], mean[i % nm]);
  
  return x;
}

