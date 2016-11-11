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
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(a))
    return NA_REAL;
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
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(a))
    return NA_REAL;
  if (sigma < 0.0 || a < 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x == R_NegInf)
    return 0.0;
  if (x == R_PosInf)
    return 1.0;
  if (sigma == 0.0)
    return R::punif(x, mu-a, mu+a, true, false);
  if (a == 0.0)
    return R::pnorm(x, mu, sigma, true, false);
  double z = x-mu;
  return sigma/(2.0*a) * (G((z+a)/sigma) - G((z-a)/sigma));
}

double rng_bhattacharjee(double mu, double sigma, double a) {
  if (ISNAN(mu) || ISNAN(sigma) || ISNAN(a))
    return NA_REAL;
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
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  dims.push_back(a.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_bhattacharjee(x[i % dims[0]], mu[i % dims[1]],
                             sigma[i % dims[2]], a[i % dims[3]]);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pbhatt(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& a,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  dims.push_back(a.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_bhattacharjee(x[i % dims[0]], mu[i % dims[1]],
                             sigma[i % dims[2]], a[i % dims[3]]);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rbhatt(
    const int& n,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& a
  ) {
  
  std::vector<int> dims;
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  dims.push_back(a.length());
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_bhattacharjee(mu[i % dims[0]], sigma[i % dims[1]], a[i % dims[2]]);
  
  return x;
}

