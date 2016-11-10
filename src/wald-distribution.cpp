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

double pdf_wald(double x, double mu, double lambda) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(lambda))
    return NA_REAL;
  if (mu <= 0.0 || lambda <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x <= 0.0 || !R_FINITE(x))
    return 0.0;
  return sqrt(lambda/(2.0*PI*pow(x, 3.0))) *
         exp((-lambda*pow(x-mu, 2.0))/(2.0*pow(mu, 2.0)*x));
}

double cdf_wald(double x, double mu, double lambda) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(lambda))
    return NA_REAL;
  if (mu <= 0.0 || lambda <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x <= 0.0)
    return 0.0;
  if (x == INFINITY)
    return 1.0;
  return Phi(sqrt(lambda/x)*(x/mu-1.0)) +
         exp((2.0*lambda)/mu) *
         Phi(-sqrt(lambda/x)*(x/mu+1.0));
}

double rng_wald(double mu, double lambda) {
  if (ISNAN(mu) || ISNAN(lambda))
    return NA_REAL;
  if (mu <= 0.0 || lambda <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
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
    bool log_prob = false
  ) {
  
  int n  = x.length();
  int nm = mu.length();
  int nl = lambda.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, nl));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_wald(x[i % n], mu[i % nm], lambda[i % nl]);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pwald(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& lambda,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = x.length();
  int nm = mu.length();
  int nl = lambda.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, nl));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_wald(x[i % n], mu[i % nm], lambda[i % nl]);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rwald(
    const int n,
    const NumericVector& mu,
    const NumericVector& lambda
  ) {
  
  int nm = mu.length();
  int nl = lambda.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_wald(mu[i % nm], lambda[i % nl]);
  
  return x;
}

