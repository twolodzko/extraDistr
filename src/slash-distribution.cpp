#include <Rcpp.h>
#include "const.h"
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
 * Location-scale slash distribution
 * 
 * Parameters:
 * mu
 * sigma > 0
 * 
 * 
 */


double pdf_slash(double x, double mu, double sigma) {
  if (sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double z = (x - mu)/sigma;
  return ((PHI_0 - phi(z))/pow(z, 2.0))/sigma;
}

double cdf_slash(double x, double mu, double sigma) {
  if (sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double z = (x - mu)/sigma;
  if (z == 0.0)
    return 0.5;
  else
    return Phi(z) - (PHI_0 - phi(z))/z;
}

double rng_slash(double mu, double sigma) {
  if (sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double z = R::norm_rand();
  double u = rng_unif();
  return z/u*sigma + mu;
}



// [[Rcpp::export]]
NumericVector cpp_dslash(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    bool log_prob = false
  ) {
  
  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_slash(x[i % n], mu[i % nm], sigma[i % ns]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pslash(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_slash(x[i % n], mu[i % nm], sigma[i % ns]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rslash(
    const int n,
    const NumericVector& mu,
    const NumericVector& sigma
  ) {
  
  int nm = mu.length();
  int ns = sigma.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_slash(mu[i % nm], sigma[i % ns]);
  
  return x;
}


