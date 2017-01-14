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
 * Wald distribution
 * 
 * Parameters:
 * mu > 0
 * lambda > 0
 * 
 * Values:
 * x > 0
 *
 * 
 */

inline double pdf_wald(double x, double mu, double lambda,
                       bool& throw_warning) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(lambda))
    return x+mu+lambda;
  if (mu <= 0.0 || lambda <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x <= 0.0 || !R_FINITE(x))
    return 0.0;
  return sqrt(lambda/(2.0*PI*pow(x, 3.0))) *
         exp((-lambda*pow(x-mu, 2.0))/(2.0*pow(mu, 2.0)*x));
}

inline double cdf_wald(double x, double mu, double lambda,
                       bool& throw_warning) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(lambda))
    return x+mu+lambda;
  if (mu <= 0.0 || lambda <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x <= 0.0)
    return 0.0;
  if (x == R_PosInf)
    return 1.0;
  return Phi(sqrt(lambda/x)*(x/mu-1.0)) +
         exp((2.0*lambda)/mu) *
         Phi(-sqrt(lambda/x)*(x/mu+1.0));
}

inline double rng_wald(double mu, double lambda, bool& throw_warning) {
  if (ISNAN(mu) || ISNAN(lambda) || mu <= 0.0 || lambda <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  double u, x, y, z;
  u = rng_unif();
  z = R::norm_rand();
  y = pow(z, 2.0);
  x = mu + (pow(mu, 2.0)*y)/(2.0*lambda) - mu/(2.0*lambda) *
      sqrt(4.0*mu*lambda*y+pow(mu, 2.0)*pow(y, 2.0));
  if (u <= mu/(mu+x))
    return x;
  else
    return pow(mu, 2.0)/x;
}


// [[Rcpp::export]]
NumericVector cpp_dwald(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& lambda,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mu.length());
  dims.push_back(lambda.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_wald(GETV(x, i), GETV(mu, i),
                    GETV(lambda, i), throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pwald(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& lambda,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mu.length());
  dims.push_back(lambda.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_wald(GETV(x, i), GETV(mu, i),
                    GETV(lambda, i), throw_warning);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rwald(
    const int& n,
    const NumericVector& mu,
    const NumericVector& lambda
  ) {
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_wald(GETV(mu, i), GETV(lambda, i),
                    throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

