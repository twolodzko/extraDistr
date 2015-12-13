#include <Rcpp.h>
#include "shared.h"
using namespace Rcpp;


/*
 * Bhattacharjee distribution
 * 
 * Parameters:
 * mu
 * sigma >= 0
 * a >= 0
 * 
 * Bhattacharjee, G.P., Pandit, S.N.N., and Mohan, R. (1963).
 * Dimensional chains involving rectangular and normal error-distributions.
 * Technometrics, 5, 404-406.
 * 
 */

double G(double x) {
  return x * Phi(x) + phi(x);
}

double pdf_bhattacharjee(double x, double mu, double sigma, double a) {
  if (sigma < 0 || a < 0)
    return NAN;
  if (sigma == 0)
    return R::dunif(x, mu-a, mu+a, false);
  if (a == 0)
    return R::dnorm(x, mu, sigma, false);
  double z = x-mu;
  return (Phi((z+a)/sigma) - Phi((z-a)/sigma)) / (2*a);
}

double cdf_bhattacharjee(double x, double mu, double sigma, double a) {
  if (sigma < 0 || a < 0)
    return NAN;
  if (sigma == 0)
    return R::punif(x, mu-a, mu+a, true, false);
  if (a == 0)
    return R::pnorm(x, mu, sigma, true, false);
  double z = x-mu;
  return sigma/(2*a) * (G((z+a)/sigma) - G((z-a)/sigma));
}

double rng_bhattacharjee(double mu, double sigma, double a) {
  if (sigma < 0 || a < 0)
    return NAN;
  if (sigma == 0)
    return R::runif(mu-a, mu+a);
  if (a == 0)
    return R::rnorm(mu, sigma);
  return R::runif(-a, +a) + R::rnorm(0, sigma) + mu;
}



// [[Rcpp::export]]
NumericVector cpp_dbhatt(NumericVector x,
                         NumericVector mu, NumericVector sigma, NumericVector a,
                         bool log_prob = false) {
  
  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int na = a.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, na));
  NumericVector p(Nmax);
  double z;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_bhattacharjee(x[i % n], mu[i % nm], sigma[i % ns], a[i % na]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pbhatt(NumericVector x,
                         NumericVector mu, NumericVector sigma, NumericVector a,
                         bool lower_tail = true, bool log_prob = false) {
  
  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int na = a.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, na));
  NumericVector p(Nmax);
  double z;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_bhattacharjee(x[i % n], mu[i % nm], sigma[i % ns], a[i % na]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rbhatt(int n,
                         NumericVector mu, NumericVector sigma, NumericVector a) {
  
  int nm = mu.length();
  int ns = sigma.length();
  int na = a.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_bhattacharjee(mu[i % nm], sigma[i % ns], a[i % na]);
  
  return x;
}

