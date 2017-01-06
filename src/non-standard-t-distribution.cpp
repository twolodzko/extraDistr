#include <Rcpp.h>

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
*  Non-standard t-distribution
*
*  Values:
*  x
*
*  Parameters:
*  nu > 0
*  mu
*  sigma > 0
*
*/

double pdf_nst(double x, double nu, double mu, double sigma) {
  if (ISNAN(x) || ISNAN(nu) || ISNAN(mu) || ISNAN(sigma))
    return NA_REAL;
  if (nu <= 0.0 || sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double z = (x - mu)/sigma;
  return R::dt(z, nu, false)/sigma;
}

double cdf_nst(double x, double nu, double mu, double sigma) {
  if (ISNAN(x) || ISNAN(nu) || ISNAN(mu) || ISNAN(sigma))
    return NA_REAL;
  if (nu <= 0.0 || sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double z = (x - mu)/sigma;
  return R::pt(z, nu, true, false);
}

double invcdf_nst(double p, double nu, double mu, double sigma) {
  if (ISNAN(p) || ISNAN(nu) || ISNAN(mu) || ISNAN(sigma))
    return NA_REAL;
  if (nu <= 0.0 || sigma <= 0.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return R::qt(p, nu, true, false)*sigma + mu;
}

double rng_nst(double nu, double mu, double sigma) {
  if (ISNAN(nu) || ISNAN(mu) || ISNAN(sigma) ||
      nu <= 0.0 || sigma <= 0.0) {
    Rcpp::warning("NAs produced");
    return NA_REAL;
  }
  return R::rt(nu)*sigma + mu;
}


// [[Rcpp::export]]
NumericVector cpp_dnst(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(nu.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_nst(x[i % dims[0]], nu[i % dims[1]], mu[i % dims[2]], sigma[i % dims[3]]);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pnst(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(nu.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_nst(x[i % dims[0]], nu[i % dims[1]], mu[i % dims[2]], sigma[i % dims[3]]);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qnst(
    const NumericVector& p,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(nu.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_nst(pp[i % dims[0]], nu[i % dims[1]], mu[i % dims[2]], sigma[i % dims[3]]);
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rnst(
    const int& n,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma
  ) {
  
  std::vector<int> dims;
  dims.push_back(nu.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_nst(nu[i % dims[0]], mu[i % dims[1]], sigma[i % dims[2]]);
  
  return x;
}

