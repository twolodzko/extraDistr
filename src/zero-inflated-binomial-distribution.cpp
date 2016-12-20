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
* Zero-inflated Poisson distribution
* 
* Parameters:
* lambda > 0
* 0 <= pi <= 1
* 
* Values:
* x >= 0
*
*/

double pdf_zib(double x, double n, double p, double pi) {
  if (ISNAN(x) || ISNAN(n) || ISNAN(p) || ISNAN(pi))
    return NA_REAL;
  if (p < 0.0 || p > 1.0 || n < 0.0 || pi < 0.0 || pi > 1.0 || !isInteger(n, false)) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0.0 || !isInteger(x) || !R_FINITE(x))
    return 0.0;
  if (x == 0.0)
    return pi + (1.0-pi) * pow(1.0-p, n);
  else
    return (1.0-pi) * R::dbinom(x, n, p, false);
}

double cdf_zib(double x, double n, double p, double pi) {
  if (ISNAN(x) || ISNAN(n) || ISNAN(p) || ISNAN(pi))
    return NA_REAL;
  if (p < 0.0 || p > 1.0 || n < 0.0 || pi < 0.0 || pi > 1.0 || !isInteger(n, false)) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  if (!R_FINITE(x))
    return 1.0;
  return pi + (1.0-pi) * R::pbinom(x, n, p, true, false);
}

double invcdf_zib(double pp, double n, double p, double pi) {
  if (ISNAN(pp) || ISNAN(n) || ISNAN(p) || ISNAN(pi))
    return NA_REAL;
  if (p < 0.0 || p > 1.0 || n < 0.0 || pi < 0.0 || pi > 1.0 ||
      !isInteger(n, false) || pp < 0.0 || pp > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (pp < pi)
    return 0.0;
  else
    return R::qbinom((pp - pi) / (1.0-pi), n, p, true, false);
}

double rng_zib(double n, double p, double pi) {
  if (ISNAN(n) || ISNAN(p) || ISNAN(pi))
    return NA_REAL;
  if (p < 0.0 || p > 1.0 || n < 0.0 || pi < 0.0 || pi > 1.0 || !isInteger(n, false)) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double u = rng_unif();
  if (u < pi)
    return 0.0;
  else
    return R::rbinom(n, p);
}


// [[Rcpp::export]]
NumericVector cpp_dzib(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(size.length());
  dims.push_back(prob.length());
  dims.push_back(pi.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_zib(x[i % dims[0]], size[i % dims[1]],
                   prob[i % dims[2]], pi[i % dims[3]]);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pzib(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(size.length());
  dims.push_back(prob.length());
  dims.push_back(pi.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_zib(x[i % dims[0]], size[i % dims[1]],
                   prob[i % dims[2]], pi[i % dims[3]]);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qzib(
    const NumericVector& p,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(size.length());
  dims.push_back(prob.length());
  dims.push_back(pi.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_zib(pp[i % dims[0]], size[i % dims[1]],
                      prob[i % dims[2]], pi[i % dims[3]]);
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rzib(
    const int& n,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi
  ) {
  
  std::vector<int> dims;
  dims.push_back(size.length());
  dims.push_back(prob.length());
  dims.push_back(pi.length());
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_zib(size[i % dims[0]], prob[i % dims[1]], pi[i % dims[2]]);
  
  return x;
}

