#include <Rcpp.h>
#include "shared.h"
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;


/*
 * Bhattacharjee distribution
 * 
 * Parameters:
 * mu
 * sigma >= 0
 * a >= 0
 * 
 * Bhattacharjee, G.P., Pandit, S.N.N., and Mohan, R. (1963).
 * Dimensional chains involving rectangular and normal error-distributions.
 * Technometrics, 5, 404-406.
 * 
 */

inline double G(double x) {
  return x * Phi(x) + phi(x);
}

inline double pdf_bhattacharjee(double x, double mu, double sigma,
                                double a, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(a))
    return x+mu+sigma+a;
#endif
  if (sigma < 0.0 || a < 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (sigma == 0.0)
    return R::dunif(x, mu-a, mu+a, false);
  if (a == 0.0)
    return R::dnorm(x, mu, sigma, false);
  double z = x-mu;
  return (Phi((z+a)/sigma) - Phi((z-a)/sigma)) / (2.0*a);
}

inline double cdf_bhattacharjee(double x, double mu, double sigma,
                                double a, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(a))
    return x+mu+sigma+a;
#endif
  if (sigma < 0.0 || a < 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x == R_NegInf)
    return 0.0;
  if (x == R_PosInf)
    return 1.0;
  if (sigma == 0.0)
    return R::punif(x, mu-a, mu+a, true, false);
  if (a == 0.0)
    return R::pnorm(x, mu, sigma, true, false);
  double z = x-mu;
  return sigma/(2.0*a) * (G((z+a)/sigma) - G((z-a)/sigma));
}

inline double rng_bhattacharjee(double mu, double sigma,
                                double a, bool& throw_warning) {
  if (ISNAN(mu) || ISNAN(sigma) || ISNAN(a) || sigma < 0.0 || a < 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  if (sigma == 0.0)
    return R::runif(mu-a, mu+a);
  if (a == 0.0)
    return R::rnorm(mu, sigma);
  return R::runif(-a, a) + R::norm_rand() * sigma + mu;
}


// [[Rcpp::export]]
NumericVector cpp_dbhatt(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& a,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), mu.length(),
                sigma.length(), a.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    mu.length(),
    sigma.length(),
    a.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_bhattacharjee(GETV(x, i), GETV(mu, i),
                             GETV(sigma, i), GETV(a, i),
                             throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pbhatt(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& a,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), mu.length(),
                sigma.length(), a.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    mu.length(),
    sigma.length(),
    a.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_bhattacharjee(GETV(x, i), GETV(mu, i),
                             GETV(sigma, i), GETV(a, i),
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
NumericVector cpp_rbhatt(
    const int& n,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& a
  ) {
  
  if (std::min({mu.length(), sigma.length(), a.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_bhattacharjee(GETV(mu, i), GETV(sigma, i),
                             GETV(a, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

