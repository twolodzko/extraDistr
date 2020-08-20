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
 *  Frechet distribution
 *
 *  Values:
 *  x > mu
 *
 *  Parameters:
 *  lambda > 0
 *  mu
 *  sigma > 0
 *
 *  z       = (x-mu)/sigma
 *  f(x)    = lambda/sigma * z^{-1-lambda} * exp(-z^-lambda)
 *  F(x)    = exp(-z^-lambda)
 *  F^-1(p) = mu + sigma * -log(p)^{-1/lambda}
 *
 */


inline double logpdf_frechet(double x, double lambda, double mu,
                             double sigma, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(mu) || ISNAN(sigma))
    return x+lambda+mu+sigma;
#endif
  if (lambda <= 0.0 || sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x <= mu)
    return R_NegInf;
  double z = (x-mu)/sigma;
  // lambda/sigma * pow(z, -1.0-lambda) * exp(-pow(z, -lambda));
  return log(lambda) - log(sigma) + log(z) * (-1.0-lambda) - exp(log(z) * -lambda);
}

inline double cdf_frechet(double x, double lambda, double mu,
                          double sigma, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(mu) || ISNAN(sigma))
    return x+lambda+mu+sigma;
#endif
  if (lambda <= 0.0 || sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x <= mu)
    return 0.0;
  double z = (x-mu)/sigma;
  return exp(-pow(z, -lambda));
}

inline double invcdf_frechet(double p, double lambda, double mu,
                             double sigma, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(lambda) || ISNAN(mu) || ISNAN(sigma))
    return p+lambda+mu+sigma;
#endif
  if (lambda <= 0.0 || sigma <= 0.0 || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  if (p == 1.0)
    return R_PosInf;
  return mu + sigma * pow(-log(p), -1.0/lambda);
}

inline double rng_frechet(double lambda, double mu,
                          double sigma, bool& throw_warning) {
  if (ISNAN(lambda) || ISNAN(mu) || ISNAN(sigma) ||
      lambda <= 0.0 || sigma <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  double u = rng_unif();
  return mu + sigma * pow(-log(u), -1.0/lambda);
}


// [[Rcpp::export]]
NumericVector cpp_dfrechet(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), lambda.length(),
                mu.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    lambda.length(),
    mu.length(),
    sigma.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_frechet(GETV(x, i), GETV(lambda, i),
                          GETV(mu, i), GETV(sigma, i),
                          throw_warning);

  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pfrechet(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), lambda.length(),
                mu.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    lambda.length(),
    mu.length(),
    sigma.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_frechet(GETV(x, i), GETV(lambda, i),
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
NumericVector cpp_qfrechet(
    const NumericVector& p,
    const NumericVector& lambda,
    const NumericVector& mu,
    const NumericVector& sigma,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), lambda.length(),
                mu.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    p.length(),
    lambda.length(),
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
    q[i] = invcdf_frechet(GETV(pp, i), GETV(lambda, i),
                          GETV(mu, i), GETV(sigma, i),
                          throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rfrechet(
    const int& n,
    const NumericVector& lambda,
    const NumericVector& mu,
    const NumericVector& sigma
  ) {
  
  if (std::min({lambda.length(), mu.length(), sigma.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }

  NumericVector x(n);
  
  bool throw_warning = false;

  for (int i = 0; i < n; i++)
    x[i] = rng_frechet(GETV(lambda, i), GETV(mu, i),
                       GETV(sigma, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

