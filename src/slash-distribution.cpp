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
using Rcpp::NumericVector;


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
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
    return NA_REAL;
  if (sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double z = (x - mu)/sigma;
  if (z == 0.0)
    return 1.0/(2.0 * SQRT_2_PI);
  return ((PHI_0 - phi(z))/pow(z, 2.0))/sigma;
}

double cdf_slash(double x, double mu, double sigma) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
    return NA_REAL;
  if (sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double z = (x - mu)/sigma;
  if (z == 0.0)
    return 0.5;
  return Phi(z) - (PHI_0 - phi(z))/z;
}

double rng_slash(double mu, double sigma) {
  if (ISNAN(mu) || ISNAN(sigma) || sigma <= 0.0) {
    Rcpp::warning("NAs produced");
    return NA_REAL;
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
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_slash(x[i % dims[0]], mu[i % dims[1]], sigma[i % dims[2]]);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pslash(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_slash(x[i % dims[0]], mu[i % dims[1]], sigma[i % dims[2]]);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rslash(
    const int& n,
    const NumericVector& mu,
    const NumericVector& sigma
  ) {
  
  std::vector<int> dims;
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_slash(mu[i % dims[0]], sigma[i % dims[1]]);
  
  return x;
}


