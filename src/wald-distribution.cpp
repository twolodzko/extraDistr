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
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(lambda))
    return x+mu+lambda;
#endif
  if (mu <= 0.0 || lambda <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x <= 0.0 || !R_FINITE(x))
    return 0.0;
  return sqrt(lambda/(2.0*M_PI*(x*x*x))) *
         exp( (-lambda*(x-mu)*(x-mu))/(2.0*(mu*mu)*x) );
}

inline double cdf_wald(double x, double mu, double lambda,
                       bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(lambda))
    return x+mu+lambda;
#endif
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
  y = z*z;
  x = mu + (mu*mu*y)/(2.0*lambda) - mu/(2.0*lambda) *
      sqrt(4.0*mu*lambda*y+(mu*mu)*(y*y));
  if (u <= mu/(mu+x))
    return x;
  else
    return (mu*mu)/x;
}


// [[Rcpp::export]]
NumericVector cpp_dwald(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& lambda,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), mu.length(), lambda.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    mu.length(),
    lambda.length()
  });
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
  
  if (std::min({x.length(), mu.length(), lambda.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    mu.length(),
    lambda.length()
  });
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
  
  if (std::min({mu.length(), lambda.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_wald(GETV(mu, i), GETV(lambda, i),
                    throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

