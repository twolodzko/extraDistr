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

using std::log1p;


/*
*  Lomax distribution
*
*  Values:
*  x > 0
*
*  Parameters:
*  lambda > 0
*  kappa > 0
*
*  f(x)    = lambda*kappa / (1+lambda*x)^(kappa+1)
*  F(x)    = 1-(1+lambda*x)^-kappa
*  F^-1(p) = ((1-p)^(-1/kappa)-1) / lambda
*
*/


inline double logpdf_lomax(double x, double lambda, double kappa,
                           bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(kappa))
    return x+lambda+kappa;
#endif
  if (lambda <= 0.0 || kappa <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x <= 0.0)
    return R_NegInf;
  // lambda*kappa / pow(1.0+lambda*x, kappa+1.0);
  return log(lambda) + log(kappa) - log1p(lambda*x)*(kappa+1.0);
}

inline double cdf_lomax(double x, double lambda, double kappa,
                        bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(kappa))
    return x+lambda+kappa;
#endif
  if (lambda <= 0.0 || kappa <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x <= 0.0)
    return 0.0;
  // 1.0 - pow(1.0+lambda*x, -kappa);
  return 1.0 - exp(log1p(lambda*x) * (-kappa));
}

inline double invcdf_lomax(double p, double lambda, double kappa,
                           bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(lambda) || ISNAN(kappa))
    return p+lambda+kappa;
#endif
  if (lambda <= 0.0 || kappa <= 0.0 || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  return (pow(1.0-p, -1.0/kappa)-1.0) / lambda;
}

inline double rng_lomax(double lambda, double kappa, bool& throw_warning) {
  if (ISNAN(lambda) || ISNAN(kappa) || lambda <= 0.0 || kappa <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  double u = rng_unif();
  return (pow(u, -1.0/kappa)-1.0) / lambda;
}


// [[Rcpp::export]]
NumericVector cpp_dlomax(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& kappa,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), lambda.length(), kappa.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    lambda.length(),
    kappa.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_lomax(GETV(x, i), GETV(lambda, i),
                        GETV(kappa, i), throw_warning);

  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_plomax(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& kappa,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), lambda.length(), kappa.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    lambda.length(),
    kappa.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_lomax(GETV(x, i), GETV(lambda, i),
                     GETV(kappa, i), throw_warning);

  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qlomax(
    const NumericVector& p,
    const NumericVector& lambda,
    const NumericVector& kappa,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), lambda.length(), kappa.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    p.length(),
    lambda.length(),
    kappa.length()
  });
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;

  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;

  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_lomax(GETV(pp, i), GETV(lambda, i),
                        GETV(kappa, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rlomax(
    const int& n,
    const NumericVector& lambda,
    const NumericVector& kappa
  ) {
  
  if (std::min({lambda.length(), kappa.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }

  NumericVector x(n);
  
  bool throw_warning = false;

  for (int i = 0; i < n; i++)
    x[i] = rng_lomax(GETV(lambda, i), GETV(kappa, i),
                     throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

