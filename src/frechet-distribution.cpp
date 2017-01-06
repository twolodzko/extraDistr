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
 *  Frechet distribution
 *
 *  Values:
 *  x > mu
 *
 *  Parameters:
 *  lambda > 0
 *  mu
 *  sigma > 0
 *
 *  z       = (x-mu)/sigma
 *  f(x)    = lambda/sigma * z^{-1-lambda} * exp(-z^-lambda)
 *  F(x)    = exp(-z^-lambda)
 *  F^-1(p) = mu + sigma * -log(p)^{-1/lambda}
 *
 */

double pdf_frechet(double x, double lambda, double mu, double sigma) {
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(mu) || ISNAN(sigma))
    return NA_REAL;
  if (lambda <= 0.0 || sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x <= mu)
    return 0.0;
  double z = (x-mu)/sigma;
  return lambda/sigma * pow(z, -1.0-lambda) * exp(-pow(z, -lambda));
}

double cdf_frechet(double x, double lambda, double mu, double sigma) {
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(mu) || ISNAN(sigma))
    return NA_REAL;
  if (lambda <= 0.0 || sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x <= mu)
    return 0.0;
  double z = (x-mu)/sigma;
  return exp(pow(-z, -lambda));
}

double invcdf_frechet(double p, double lambda, double mu, double sigma) {
  if (ISNAN(p) || ISNAN(lambda) || ISNAN(mu) || ISNAN(sigma))
    return NA_REAL;
  if (lambda <= 0.0 || sigma <= 0.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (p == 1.0)
    return R_PosInf;
  return mu + sigma * pow(-log(p), -1.0/lambda);
}

double rng_frechet(double lambda, double mu, double sigma) {
  if (ISNAN(lambda) || ISNAN(mu) || ISNAN(sigma) ||
      lambda <= 0.0 || sigma <= 0.0) {
    Rcpp::warning("NAs produced");
    return NA_REAL;
  }
  double u = rng_unif();
  return mu + sigma * pow(-log(u), -1.0/lambda);
}


// [[Rcpp::export]]
NumericVector cpp_dfrechet(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(lambda.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_frechet(x[i % dims[0]], lambda[i % dims[1]],
                       mu[i % dims[2]], sigma[i % dims[3]]);

  if (log_prob)
    p = Rcpp::log(p);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pfrechet(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(lambda.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_frechet(x[i % dims[0]], lambda[i % dims[1]],
                       mu[i % dims[2]], sigma[i % dims[3]]);

  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qfrechet(
    const NumericVector& p,
    const NumericVector& lambda,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(lambda.length());
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
    q[i] = invcdf_frechet(pp[i % dims[0]], lambda[i % dims[1]],
                          mu[i % dims[2]], sigma[i % dims[3]]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rfrechet(
    const int& n,
    const NumericVector& lambda,
    const NumericVector& mu,
    const NumericVector& sigma
  ) {

  std::vector<int> dims;
  dims.push_back(lambda.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_frechet(lambda[i % dims[0]], mu[i % dims[1]], sigma[i % dims[2]]);

  return x;
}

