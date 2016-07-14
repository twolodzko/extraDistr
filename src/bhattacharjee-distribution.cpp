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
  if (sigma < 0.0 || a < 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (sigma == 0.0)
    return R::dunif(x, mu-a, mu+a, false);
  if (a == 0.0)
    return R::dnorm(x, mu, sigma, false);
  double z = x-mu;
  return (Phi((z+a)/sigma) - Phi((z-a)/sigma)) / (2.0*a);
}

double cdf_bhattacharjee(double x, double mu, double sigma, double a) {
  if (sigma < 0.0 || a < 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x == -INFINITY)
    return 0.0;
  if (x == INFINITY)
    return 1.0;
  if (sigma == 0.0)
    return R::punif(x, mu-a, mu+a, true, false);
  if (a == 0.0)
    return R::pnorm(x, mu, sigma, true, false);
  double z = x-mu;
  return sigma/(2.0*a) * (G((z+a)/sigma) - G((z-a)/sigma));
}

double rng_bhattacharjee(double mu, double sigma, double a) {
  if (sigma < 0.0 || a < 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (sigma == 0.0)
    return R::runif(mu-a, mu+a);
  if (a == 0.0)
    return R::rnorm(mu, sigma);
  return R::runif(-a, a) + R::norm_rand() * sigma + mu;
}



// [[Rcpp::export]]
NumericVector cpp_dbhatt(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& a,
    bool log_prob = false
  ) {
  
  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int na = a.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, na));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_bhattacharjee(x[i % n], mu[i % nm], sigma[i % ns], a[i % na]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pbhatt(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& a,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int na = a.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, na));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_bhattacharjee(x[i % n], mu[i % nm], sigma[i % ns], a[i % na]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rbhatt(
    const int n,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& a
  ) {
  
  int nm = mu.length();
  int ns = sigma.length();
  int na = a.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_bhattacharjee(mu[i % nm], sigma[i % ns], a[i % na]);
  
  return x;
}

