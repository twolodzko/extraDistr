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
 * Discrete uniform distribution
 * 
 * Values:
 * a <= x <= b
 * 
 * f(x) = 1/(b-a+1)
 * F(x) = (floor(x)-a+1)/b-a+1
 *  
 */


double pmf_dunif(double x, double min, double max) {
  if (ISNAN(x) || ISNAN(min) || ISNAN(max))
    return NA_REAL;
  if (min > max || !R_FINITE(min) || !R_FINITE(max) ||
      !isInteger(min, false) || !isInteger(max, false)) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < min || x > max || !isInteger(x))
    return 0.0;
  return 1.0/(max-min+1.0);
}


double cdf_dunif(double x, double min, double max) {
  if (ISNAN(x) || ISNAN(min) || ISNAN(max))
    return NA_REAL;
  if (min > max || !R_FINITE(min) || !R_FINITE(max) ||
      !isInteger(min, false) || !isInteger(max, false)) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < min)
    return 0.0;
  else if (x >= max)
    return 1.0;
  return (floor(x)-min+1.0)/(max-min+1.0);
}

double invcdf_dunif(double p, double min, double max) {
  if (ISNAN(p) || ISNAN(min) || ISNAN(max))
    return NA_REAL;
  if (min > max || !R_FINITE(min) || !R_FINITE(max) ||
      !isInteger(min, false) || !isInteger(max, false) ||
      p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (p == 0 || min == max)
    return min;
  return ceil( p*(max-min+1.0)+min-1.0 );
}

double rng_dunif(double min, double max) {
  if (ISNAN(min) || ISNAN(max) ||
      min > max || !R_FINITE(min) || !R_FINITE(max) ||
      !isInteger(min, false) || !isInteger(max, false)) {
    Rcpp::warning("NAs produced");
    return NA_REAL;
  }
  if (min == max)
    return min;
  return ceil(R::runif(min - 1.0, max));
}


// [[Rcpp::export]]
NumericVector cpp_ddunif(
    const NumericVector& x,
    const NumericVector& min,
    const NumericVector& max,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(min.length());
  dims.push_back(max.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pmf_dunif(x[i % dims[0]], min[i % dims[1]], max[i % dims[2]]);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pdunif(
    const NumericVector& x,
    const NumericVector& min,
    const NumericVector& max,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(min.length());
  dims.push_back(max.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_dunif(x[i % dims[0]], min[i % dims[1]], max[i % dims[2]]);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qdunif(
    const NumericVector& p,
    const NumericVector& min,
    const NumericVector& max,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(min.length());
  dims.push_back(max.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_dunif(pp[i % dims[0]], min[i % dims[1]], max[i % dims[2]]);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rdunif(
    const int& n,
    const NumericVector& min,
    const NumericVector& max
  ) {
  
  std::vector<int> dims;
  dims.push_back(min.length());
  dims.push_back(max.length());
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_dunif(min[i % dims[0]], max[i % dims[1]]);
  
  return x;
}

