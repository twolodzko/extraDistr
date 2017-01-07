#include <Rcpp.h>
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
 *  Laplace distribution
 *
 *  Values:
 *  x
 *
 *  Parameters:
 *  mu
 *  sigma > 0
 *
 *  z = (x-mu)/sigma
 *  f(x)    = 1/(2*sigma) * exp(-|z|)
 *  F(x)    = { 1/2 * exp(z)                 if   x < mu
 *            { 1 - 1/2 * exp(z)             otherwise
 *  F^-1(p) = { mu + sigma * log(2*p)        if p <= 0.5
 *            { mu + sigma * log(2*(1-p))    otherwise
 *
 */

double pdf_laplace(double x, double mu, double sigma) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
    return NA_REAL;
  if (sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double z = abs(x-mu)/sigma;
  return exp(-z)/(2.0*sigma);
}

double cdf_laplace(double x, double mu, double sigma) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
    return NA_REAL;
  if (sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double z = (x-mu)/sigma;
  if (x < mu)
    return exp(z)/2.0;
  else
    return 1.0 - exp(-z)/2.0;
}

double invcdf_laplace(double p, double mu, double sigma) {
  if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma))
    return NA_REAL;
  if (sigma <= 0.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (p < 0.5)
    return mu + sigma * log(2.0*p);
  else
    return mu - sigma * log(2.0*(1.0-p));
}

double rng_laplace(double mu, double sigma) {
  if (ISNAN(mu) || ISNAN(sigma) || sigma <= 0.0) {
    Rcpp::warning("NAs produced");
    return NA_REAL;
  }
  // this is slower
  // double u = R::runif(-0.5, 0.5);
  // return mu + sigma * R::sign(u) * log(1.0 - 2.0*abs(u));
  double u = R::exp_rand();
  double s = rng_sign();
  return u*s * sigma + mu;
}


// [[Rcpp::export]]
NumericVector cpp_dlaplace(
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
    p[i] = pdf_laplace(x[i % dims[0]], mu[i % dims[1]], sigma[i % dims[2]]);

  if (log_prob)
    p = Rcpp::log(p);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_plaplace(
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
    p[i] = cdf_laplace(x[i % dims[0]], mu[i % dims[1]], sigma[i % dims[2]]);

  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qlaplace(
    const NumericVector& p,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);

  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_laplace(pp[i % dims[0]], mu[i % dims[1]], sigma[i % dims[2]]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rlaplace(
    const int& n,
    const NumericVector& mu,
    const NumericVector& sigma
  ) {

  std::vector<int> dims;
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_laplace(mu[i % dims[0]], sigma[i % dims[1]]);

  return x;
}

