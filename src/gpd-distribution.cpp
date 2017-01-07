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
*  Generalized Pareto distribution
*
*  Values:
*  x
*
*  Parameters:
*  mu
*  sigma > 0
*  xi
*
*  z = (x-mu)/sigma
*  where 1+xi*z > 0
*
*  f(x)    = { (1+xi*z)^{-(xi+1)/xi}/sigma       if xi != 0
*            { exp(-z)/sigma                     otherwise
*  F(x)    = { 1-(1+xi*z)^{-1/xi}                if xi != 0
*            { 1-exp(-z)                         otherwise
*  F^-1(p) = { mu + sigma * ((1-p)^{-xi}-1)/xi   if xi != 0
*            { mu - sigma * log(1-p)             otherwise
*
*/

double pdf_gpd(double x, double mu, double sigma, double xi) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(xi))
    return NA_REAL;
  if (sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double z = (x-mu)/sigma;
  if (xi != 0.0) {
    if (x >= mu)
      return pow(1.0+xi*z, -(xi+1.0)/xi)/sigma;
    else
      return 0.0;
  } else {
    if (x >= mu && x <= (mu - sigma/xi))
      return exp(-z)/sigma;
    else
      return 0.0;
  }
}

double cdf_gpd(double x, double mu, double sigma, double xi) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(xi))
    return NA_REAL;
  if (sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double z = (x-mu)/sigma;
  if (xi != 0.0) {
    if (x >= mu)
      return 1.0 - pow(1.0+xi*z, -1.0/xi);
    else
      return 0.0;
  } else {
    if (x >= mu && x <= (mu - sigma/xi))
      return 1.0 - exp(-z);
    else
      return 0.0;
  }
}

double invcdf_gpd(double p, double mu, double sigma, double xi) {
  if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma) || ISNAN(xi))
    return NA_REAL;
  if (sigma <= 0.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (xi != 0.0)
    return mu + sigma * (pow(1.0-p, -xi)-1.0)/xi;
  else
    return mu - sigma * log(1.0-p);
}

double rng_gpd(double mu, double sigma, double xi) {
  if (ISNAN(mu) || ISNAN(sigma) || ISNAN(xi) || sigma <= 0.0) {
    Rcpp::warning("NAs produced");
    return NA_REAL;
  }
  double u = rng_unif();
  if (xi != 0.0)
    return mu + sigma * (pow(u, -xi)-1.0)/xi;
  else
    return mu - sigma * log(u);
}


// [[Rcpp::export]]
NumericVector cpp_dgpd(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& xi,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  dims.push_back(xi.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_gpd(x[i % dims[0]], mu[i % dims[1]],
                   sigma[i % dims[2]], xi[i % dims[3]]);

  if (log_prob)
    p = Rcpp::log(p);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pgpd(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& xi,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  dims.push_back(xi.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_gpd(x[i % dims[0]], mu[i % dims[1]],
                   sigma[i % dims[2]], xi[i % dims[3]]);

  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qgpd(
    const NumericVector& p,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& xi,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  dims.push_back(xi.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);

  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_gpd(pp[i % dims[0]], mu[i % dims[1]],
                      sigma[i % dims[2]], xi[i % dims[3]]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rgpd(
    const int& n,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& xi
  ) {

  std::vector<int> dims;
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  dims.push_back(xi.length());
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_gpd(mu[i % dims[0]], sigma[i % dims[1]], xi[i % dims[2]]);

  return x;
}

