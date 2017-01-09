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

double pdf_zinb(double x, double r, double p,
                double pi, bool& throw_warning) {
  if (ISNAN(x) || ISNAN(r) || ISNAN(p) || ISNAN(pi))
    return x+r+p+pi;
  if (p < 0.0 || p > 1.0 || r < 0.0 || pi < 0.0 || pi > 1.0 ||
      !isInteger(r, false)) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0 || !isInteger(x) || !R_FINITE(x))
    return 0.0;
  if (x == 0.0)
    return pi + (1.0-pi) * pow(p, r);
  else
    return (1.0-pi) * R::dnbinom(x, r, p, false);
}

double cdf_zinb(double x, double r, double p,
                double pi, bool& throw_warning) {
  if (ISNAN(x) || ISNAN(r) || ISNAN(p) || ISNAN(pi))
    return x+r+p+pi;
  if (p < 0.0 || p > 1.0 || r < 0.0 || pi < 0.0 || pi > 1.0 ||
      !isInteger(r, false)) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  if (!R_FINITE(x))
    return 1.0;
  return pi + (1.0-pi) * R::pnbinom(x, r, p, true, false);
}

double invcdf_zinb(double pp, double r, double p,
                   double pi, bool& throw_warning) {
  if (ISNAN(pp) || ISNAN(r) || ISNAN(p) || ISNAN(pi))
    return pp+r+p+pi;
  if (p < 0.0 || p > 1.0 || r < 0.0 || pi < 0.0 || pi > 1.0 ||
      !isInteger(r, false) || pp < 0.0 || pp > 1.0) {
    throw_warning = true;
    return NAN;
  }
  if (pp < pi)
    return 0.0;
  else
    return R::qnbinom((pp - pi) / (1.0-pi), r, p, true, false);
}

double rng_zinb(double r, double p, double pi,
                bool& throw_warning) {
  if (ISNAN(r) || ISNAN(p) || ISNAN(pi) ||
      p < 0.0 || p > 1.0 || r < 0.0 ||
      pi < 0.0 || pi > 1.0 || !isInteger(r, false)) {
    throw_warning = true;
    return NA_REAL;
  }
  double u = rng_unif();
  if (u < pi)
    return 0.0;
  else
    return R::rnbinom(r, p);
}


// [[Rcpp::export]]
NumericVector cpp_dzinb(
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
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_zinb(x[i % dims[0]], size[i % dims[1]],
                    prob[i % dims[2]], pi[i % dims[3]],
                    throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pzinb(
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
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_zinb(x[i % dims[0]], size[i % dims[1]],
                    prob[i % dims[2]], pi[i % dims[3]],
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
NumericVector cpp_qzinb(
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
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_zinb(pp[i % dims[0]], size[i % dims[1]],
                       prob[i % dims[2]], pi[i % dims[3]],
                       throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rzinb(
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
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_zinb(size[i % dims[1]], prob[i % dims[1]],
                    pi[i % dims[2]], throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

