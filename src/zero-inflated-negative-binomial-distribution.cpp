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

double pdf_zinb(double x, double r, double p, double pi) {
  if (ISNAN(x) || ISNAN(r) || ISNAN(p) || ISNAN(pi))
    return NAN;
  // if (p < 0.0 || p > 1.0 || r < 0.0 || pi < 0.0 || pi > 1.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }
  if (x < 0.0 || !isInteger(x) || std::isinf(x))
    return 0.0;
  if (x == 0.0)
    return pi + (1.0-pi) * pow(p, r);
  else
    return (1.0-pi) * R::dnbinom(x, r, p, false);
}

double cdf_zinb(double x, double r, double p, double pi) {
  if (ISNAN(x) || ISNAN(r) || ISNAN(p) || ISNAN(pi))
    return NAN;
  // if (p < 0.0 || p > 1.0 || r < 0.0 || pi < 0.0 || pi > 1.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }
  if (x < 0.0)
    return 0.0;
  if (std::isinf(x))
    return 1.0;
  return pi + (1.0-pi) * R::pnbinom(x, r, p, true, false);
}

double invcdf_zinb(double pp, double r, double p, double pi) {
  if (ISNAN(pp) || ISNAN(r) || ISNAN(p) || ISNAN(pi))
    return NAN;
  // if (p < 0.0 || p > 1.0 || r < 0.0 || pi < 0.0 || pi > 1.0 || pp < 0.0 || pp > 1.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }
  if (pp < pi)
    return 0.0;
  else
    return R::qnbinom((pp - pi) / (1.0-pi), r, p, true, false);
}

double rng_zinb(double r, double p, double pi) {
  if (ISNAN(r) || ISNAN(p) || ISNAN(pi))
    return NAN;
  // if (p < 0.0 || p > 1.0 || r < 0.0 || pi < 0.0 || pi > 1.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }
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
    bool log_prob = false
  ) {
  
  int n  = x.length();
  int npi = pi.length();
  int ns = size.length();
  int np = prob.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, npi, ns, np));
  NumericVector p(Nmax);
  NumericVector size_n = nonneg_or_nan(size);
  NumericVector prob_n = zeroone_or_nan(prob);
  NumericVector pi_n = zeroone_or_nan(pi);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_zinb(x[i % n], size_n[i % ns], prob_n[i % np], pi_n[i % npi]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pzinb(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = x.length();
  int npi = pi.length();
  int ns = size.length();
  int np = prob.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, npi, ns, np));
  NumericVector p(Nmax);
  NumericVector size_n = nonneg_or_nan(size);
  NumericVector prob_n = zeroone_or_nan(prob);
  NumericVector pi_n = zeroone_or_nan(pi);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_zinb(x[i % n], size_n[i % ns], prob_n[i % np], pi_n[i % np]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qzinb(
    const NumericVector& p,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = p.length();
  int npi = pi.length();
  int ns = size.length();
  int np = prob.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, npi, ns, np));
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  NumericVector size_n = nonneg_or_nan(size);
  NumericVector prob_n = zeroone_or_nan(prob);
  NumericVector pi_n = zeroone_or_nan(pi);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1.0 - pp[i];
  
  pp = zeroone_or_nan(pp);
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_zinb(pp[i % n], size_n[i % ns], prob_n[i % np], pi_n[i % np]);
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rzinb(
    const int n,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& pi
  ) {
  
  int npi = pi.length();
  int ns = size.length();
  int np = prob.length();
  NumericVector x(n);
  NumericVector size_n = nonneg_or_nan(size);
  NumericVector prob_n = zeroone_or_nan(prob);
  NumericVector pi_n = zeroone_or_nan(pi);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_zinb(size_n[i % ns], prob_n[i % np], pi_n[i % npi]);
  
  return x;
}

