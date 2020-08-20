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
*  Non-standard t-distribution
*
*  Values:
*  x
*
*  Parameters:
*  nu > 0
*  mu
*  sigma > 0
*
*/

inline double pdf_lst(double x, double nu, double mu, double sigma,
                      bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(nu) || ISNAN(mu) || ISNAN(sigma))
    return x+nu+mu+sigma;
#endif
  if (nu <= 0.0 || sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  double z = (x - mu)/sigma;
  return R::dt(z, nu, false)/sigma;
}

inline double cdf_lst(double x, double nu, double mu, double sigma,
                      bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(nu) || ISNAN(mu) || ISNAN(sigma))
    return x+nu+mu+sigma;
#endif
  if (nu <= 0.0 || sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  double z = (x - mu)/sigma;
  return R::pt(z, nu, true, false);
}

inline double invcdf_lst(double p, double nu, double mu, double sigma,
                         bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(nu) || ISNAN(mu) || ISNAN(sigma))
    return p+nu+mu+sigma;
#endif
  if (nu <= 0.0 || sigma <= 0.0 || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  return R::qt(p, nu, true, false)*sigma + mu;
}

inline double rng_lst(double nu, double mu, double sigma,
                      bool& throw_warning) {
  if (ISNAN(nu) || ISNAN(mu) || ISNAN(sigma) ||
      nu <= 0.0 || sigma <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  return R::rt(nu)*sigma + mu;
}


// [[Rcpp::export]]
NumericVector cpp_dlst(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), nu.length(),
                mu.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    nu.length(),
    mu.length(),
    sigma.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_lst(GETV(x, i), GETV(nu, i),
                   GETV(mu, i), GETV(sigma, i),
                   throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_plst(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), nu.length(),
                mu.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    nu.length(),
    mu.length(),
    sigma.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_lst(GETV(x, i), GETV(nu, i),
                   GETV(mu, i), GETV(sigma, i),
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
NumericVector cpp_qlst(
    const NumericVector& p,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), nu.length(),
                mu.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    p.length(),
    nu.length(),
    mu.length(),
    sigma.length()
  });
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_lst(GETV(pp, i), GETV(nu, i),
                      GETV(mu, i), GETV(sigma, i),
                      throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rlst(
    const int& n,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma
  ) {
  
  if (std::min({nu.length(), mu.length(), sigma.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_lst(GETV(nu, i), GETV(mu, i),
                   GETV(sigma, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

