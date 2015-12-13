#include <Rcpp.h>
#include "shared.h"
using namespace Rcpp;

/*
 * Location-scale slash distribution
 * 
 * Parameters:
 * mu
 * sigma > 0
 * 
 * 
 */


double pdf_slash(double x, double mu, double sigma) {
  if (sigma <= 0)
    return NAN;
  double z = (x - mu)/sigma;
  return ((phi(0) - phi(z))/std::pow(z, 2))/sigma;
}

double cdf_slash(double x, double mu, double sigma) {
  if (sigma <= 0)
    return NAN;
  double z = (x - mu)/sigma;
  if (z == 0)
    return 0.5;
  else
    return Phi(z) - (phi(0) - phi(z))/z;
}

double rng_slash(double mu, double sigma) {
  if (sigma <= 0)
    return NAN;
  double z = R::rnorm(0, 1);
  double u = R::runif(0, 1);
  return z/u*sigma + mu;
}



// [[Rcpp::export]]
NumericVector cpp_dslash(NumericVector x,
                       NumericVector mu, NumericVector sigma,
                       bool log_prob = false) {
  
  double z;
  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_slash(x[i % n], mu[i % nm], sigma[i % ns]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pslash(NumericVector x,
                       NumericVector mu, NumericVector sigma,
                       bool lower_tail = true, bool log_prob = false) {
  
  double z;
  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_slash(x[i % n], mu[i % nm], sigma[i % ns]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rslash(int n,
                         NumericVector mu, NumericVector sigma) {
  
  double u;
  int nm = mu.length();
  int ns = sigma.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_slash(mu[i % nm], sigma[i % ns]);
  
  return x;
}


