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
 *  Gumbel distribution
 *
 *  Values:
 *  x
 *
 *  Parameters:
 *  mu
 *  sigma > 0
 *
 *  z       = (x-mu)/sigma
 *  f(x)    = 1/sigma * exp(-(z+exp(-z)))
 *  F(x)    = exp(-exp(-z))
 *  F^-1(p) = mu - sigma * log(-log(p))
 *
 */

inline double logpdf_gumbel(double x, double mu, double sigma,
                            bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
    return x+mu+sigma;
#endif
  if (sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (!R_FINITE(x))
    return R_NegInf;
  double z = (x-mu)/sigma;
  // exp(-(z+exp(-z)))/sigma;
  return -(z+exp(-z)) - log(sigma);
}


inline double cdf_gumbel(double x, double mu, double sigma,
                         bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
    return x+mu+sigma;
#endif
  if (sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  double z = (x-mu)/sigma;
  return exp(-exp(-z));
}

inline double invcdf_gumbel(double p, double mu, double sigma,
                            bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma))
    return p+mu+sigma;
#endif
  if (sigma <= 0.0 || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  return mu - sigma * log(-log(p));
}

inline double rng_gumbel(double mu, double sigma,
                         bool& throw_warning) {
  if (ISNAN(mu) || ISNAN(sigma) || sigma <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  double u = R::exp_rand(); // -log(rng_unif())
  return mu - sigma * log(u);
}


// [[Rcpp::export]]
NumericVector cpp_dgumbel(
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
    p[i] = logpdf_gumbel(GETV(x, i), GETV(mu, i),
                         GETV(sigma, i), throw_warning);

  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pgumbel(
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
    p[i] = cdf_gumbel(GETV(x, i), GETV(mu, i),
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
NumericVector cpp_qgumbel(
    const NumericVector& p,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), mu.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    p.length(),
    mu.length(),
    sigma.length()
  });
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;

  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_gumbel(GETV(pp, i), GETV(mu, i),
                         GETV(sigma, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rgumbel(
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
    x[i] = rng_gumbel(GETV(mu, i), GETV(sigma, i),
                      throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

