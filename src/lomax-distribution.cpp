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
*  Lomax distribution
*
*  Values:
*  x > 0
*
*  Parameters:
*  lambda > 0
*  kappa > 0
*
*  f(x)    = lambda*kappa / (1+lambda*x)^(kappa+1)
*  F(x)    = 1-(1+lambda*x)^-kappa
*  F^-1(p) = ((1-p)^(-1/kappa)-1) / lambda
*
*/

double pdf_lomax(double x, double lambda, double kappa) {
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(kappa))
    return NA_REAL;
  if (lambda <= 0.0 || kappa <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x <= 0.0)
    return 0.0;
  return lambda*kappa / pow(1.0+lambda*x, kappa+1.0);
}

double logpdf_lomax(double x, double lambda, double kappa) {
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(kappa))
    return NA_REAL;
  if (lambda <= 0.0 || kappa <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x <= 0.0)
    return R_NegInf;
  return log(lambda) + log(kappa) - log(1.0+lambda*x)*(kappa+1.0);
}

double cdf_lomax(double x, double lambda, double kappa) {
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(kappa))
    return NA_REAL;
  if (lambda <= 0.0 || kappa <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x <= 0.0)
    return 0.0;
  return 1.0 - pow(1.0+lambda*x, -kappa);
}

double invcdf_lomax(double p, double lambda, double kappa) {
  if (ISNAN(p) || ISNAN(lambda) || ISNAN(kappa))
    return NA_REAL;
  if (lambda <= 0.0 || kappa <= 0.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return (pow(1.0-p, -1.0/kappa)-1.0) / lambda;
}


// [[Rcpp::export]]
NumericVector cpp_dlomax(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& kappa,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(lambda.length());
  dims.push_back(kappa.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_lomax(x[i % dims[0]], lambda[i % dims[1]], kappa[i % dims[2]]);

  if (!log_prob)
    p = Rcpp::exp(p);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_plomax(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& kappa,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(lambda.length());
  dims.push_back(kappa.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_lomax(x[i % dims[0]], lambda[i % dims[1]], kappa[i % dims[2]]);

  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qlomax(
    const NumericVector& p,
    const NumericVector& lambda,
    const NumericVector& kappa,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(lambda.length());
  dims.push_back(kappa.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);

  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;

  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_lomax(pp[i % dims[0]], lambda[i % dims[1]], kappa[i % dims[2]]);

  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rlomax(
    const int& n,
    const NumericVector& lambda,
    const NumericVector& kappa
  ) {

  double u;
  std::vector<int> dims;
  dims.push_back(lambda.length());
  dims.push_back(kappa.length());
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    u = rng_unif();
    x[i] = invcdf_lomax(u, lambda[i % dims[0]], kappa[i % dims[1]]);
  }

  return x;
}

