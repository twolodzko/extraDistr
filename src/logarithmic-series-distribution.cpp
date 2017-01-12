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
*  Logarithmic Series distribution
*
*  Values:
*  x
*
*  Parameters:
*  0 < theta < 1
*
*  f(x) = (-1/log(1-theta)*theta^x) / x
*  F(x) = -1/log(1-theta) * sum((theta^x)/x)
*
*/


double pdf_lgser(double x, double theta, bool& throw_warnin) {
  if (ISNAN(x) || ISNAN(theta))
    return x+theta;
  if (theta <= 0.0 || theta >= 1.0) {
    throw_warnin = true;
    return NAN;
  }
  if (!isInteger(x) || x < 1.0)
    return 0.0;
  double a = -1.0/log(1.0 - theta);
  return a * pow(theta, x) / x;
}


double cdf_lgser(double x, double theta, bool& throw_warnin) {
  if (ISNAN(x) || ISNAN(theta))
    return x+theta;
  if (theta <= 0.0 || theta >= 1.0) {
    throw_warnin = true;
    return NAN;
  }
  if (x < 1.0)
    return 0.0;
  if (!R_FINITE(x))
    return 1.0;
  
  double a = -1.0/log(1.0 - theta);
  double b = 0.0;
  
  for (double k = 1.0; k <= x; k += 1.0)
    b += pow(theta, k) / k;
  
  return a * b;
}

double invcdf_lgser(double p, double theta, bool& throw_warnin) {
  if (ISNAN(p) || ISNAN(theta))
    return p+theta;
  if (theta <= 0.0 || theta >= 1.0 || !VALID_PROB(p)) {
    throw_warnin = true;
    return NAN;
  }
  if (p == 0.0)
    return 1.0;
  if (p == 1.0)
    return R_PosInf;
  
  double pk = -theta/log(1.0 - theta);
  double k = 1.0;
  
  while (p > pk) {
    p -= pk;
    pk *= theta * k/(k+1.0);
    k += 1.0;
  }
  
  return k;
}

double rng_lgser(double theta, bool& throw_warnin) {
  if (ISNAN(theta) || theta <= 0.0 || theta >= 1.0) {
    throw_warnin = true;
    return NA_REAL;
  }

  double u = rng_unif();
  double pk = -theta/log(1.0 - theta);
  double k = 1.0;
  
  while (u > pk) {
    u -= pk;
    pk *= theta * k/(k+1.0);
    k += 1.0;
  }
  
  return k;
}


// [[Rcpp::export]]
NumericVector cpp_dlgser(
    const NumericVector& x,
    const NumericVector& theta,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(theta.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_lgser(x[i % dims[0]], theta[i % dims[1]],
                     throw_warning);
 
 if (log_prob)
   p = Rcpp::log(p);
 
 if (throw_warning)
   Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_plgser(
    const NumericVector& x,
    const NumericVector& theta,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(theta.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_lgser(x[i % dims[0]], theta[i % dims[1]],
                     throw_warning);

  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);

  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qlgser(
    const NumericVector& p,
    const NumericVector& theta,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(theta.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_lgser(pp[i % dims[0]], theta[i % dims[1]],
                        throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rlgser(
    const int& n,
    const NumericVector& theta
  ) {

  int dims = theta.length();
  NumericVector x(n);
  
  bool throw_warning = false;

  for (int i = 0; i < n; i++)
    x[i] = rng_lgser(theta[i % dims], throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

