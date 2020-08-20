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
 * Location-scale slash distribution
 * 
 * Parameters:
 * mu
 * sigma > 0
 * 
 * 
 */


inline double pdf_slash(double x, double mu, double sigma,
                        bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
    return x+mu+sigma;
#endif
  if (sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  double z = (x - mu)/sigma;
  if (z == 0.0)
    return 0.19947114020071635; // 1.0/(2.0 * SQRT_2_PI);
  return ((PHI_0 - phi(z))/(z*z))/sigma;
}

inline double cdf_slash(double x, double mu, double sigma,
                        bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
    return x+mu+sigma;
#endif
  if (sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  double z = (x - mu)/sigma;
  if (z == 0.0)
    return 0.5;
  return Phi(z) - (PHI_0 - phi(z))/z;
}

inline double rng_slash(double mu, double sigma,
                        bool& throw_warning) {
  if (ISNAN(mu) || ISNAN(sigma) || sigma <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  double z = R::norm_rand();
  double u = rng_unif();
  return z/u*sigma + mu;
}


// [[Rcpp::export]]
NumericVector cpp_dslash(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), mu.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    mu.length(),
    sigma.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_slash(GETV(x, i), GETV(mu, i),
                     GETV(sigma, i), throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pslash(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), mu.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    mu.length(),
    sigma.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_slash(GETV(x, i), GETV(mu, i),
                     GETV(sigma, i), throw_warning);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rslash(
    const int& n,
    const NumericVector& mu,
    const NumericVector& sigma
  ) {
  
  if (std::min({mu.length(), sigma.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_slash(GETV(mu, i), GETV(sigma, i),
                     throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}


