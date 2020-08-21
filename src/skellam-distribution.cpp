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
 * Skellam distribution
 * 
 * mu1 >= 0
 * mu2 >= 0
 * 
 */

inline double pmf_skellam(double x, double mu1, double mu2,
                          bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu1) || ISNAN(mu2))
    return x+mu1+mu2;
#endif
  if (mu1 < 0.0 || mu2 < 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (!isInteger(x) || !R_FINITE(x))
    return 0.0;
  return exp(-(mu1+mu2)) * pow(mu1/mu2, x/2.0) *
    R::bessel_i(2.0*sqrt(mu1*mu2), x, 1.0);
}

inline double rng_skellam(double mu1, double mu2,
                          bool& throw_warning) {
  if (ISNAN(mu1) || ISNAN(mu2) || mu1 < 0.0 || mu2 < 0.0) {
    throw_warning = true;
    return NA_REAL;
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
  
  if (std::min({x.length(), mu1.length(), mu2.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    mu1.length(),
    mu2.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pmf_skellam(GETV(x, i), GETV(mu1, i),
                       GETV(mu2, i), throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rskellam(
    const int& n,
    const NumericVector& mu1,
    const NumericVector& mu2
  ) {
  
  if (std::min({mu1.length(), mu2.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_skellam(GETV(mu1, i), GETV(mu2, i),
                       throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

