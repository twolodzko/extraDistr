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
 * Skellam distribution
 * 
 * mu1 >= 0
 * mu2 >= 0
 * 
 */

double pmf_skellam(double x, double mu1, double mu2) {
  if (ISNAN(x) || ISNAN(mu1) || ISNAN(mu2))
    return NA_REAL;
  if (mu1 < 0.0 || mu2 < 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(x) || !R_finite(x))
    return 0.0;
  return exp(-(mu1+mu2)) * pow(mu1/mu2, x/2.0) * R::bessel_i(2.0*sqrt(mu1*mu2), x, 1.0);
}

double rng_skellam(double mu1, double mu2) {
  if (ISNAN(mu1) || ISNAN(mu2))
    return NA_REAL;
  if (mu1 < 0.0 || mu2 < 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return R::rpois(mu1) - R::rpois(mu2);
}



// [[Rcpp::export]]
NumericVector cpp_dskellam(
    const NumericVector& x,
    const NumericVector& mu1,
    const NumericVector& mu2,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mu1.length());
  dims.push_back(mu2.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pmf_skellam(x[i % dims[0]], mu1[i % dims[1]], mu2[i % dims[2]]);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rskellam(
    const int& n,
    const NumericVector& mu1,
    const NumericVector& mu2
  ) {
  
  std::vector<int> dims;
  dims.push_back(mu1.length());
  dims.push_back(mu2.length());
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_skellam(mu1[i % dims[0]], mu2[i % dims[1]]);
  
  return x;
}

