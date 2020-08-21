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
*  Inverse-Gamma distribution
*
*  Values:
*  x
*
*  Parameters:
*  alpha > 0
*  beta > 0
*
*  f(k) = (x^(-alpha-1) * exp(-1/(beta*x))) / (Gamma(alpha)*beta^alpha)
*  F(x) = gamma(alpha, 1/(beta*x)) / Gamma(alpha)
*
*  V. Witkovsky (2001) Computing the distribution of a linear
*  combination of inverted gamma variables, Kybernetika 37(1), 79-90
*
*/


inline double logpdf_invgamma(double x, double alpha, double beta,
                              bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(alpha) || ISNAN(beta))
    return x+alpha+beta;
#endif
  if (alpha <= 0.0 || beta <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x <= 0.0)
    return R_NegInf;
  
  return log(beta) * -alpha - R::lgammafn(alpha) + log(x) *
           (-alpha-1.0) - 1.0/(beta*x);
}

inline double cdf_invgamma(double x, double alpha, double beta,
                           bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(alpha) || ISNAN(beta))
    return x+alpha+beta;
#endif
  if (alpha <= 0.0 || beta <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x <= 0.0)
    return 0.0;
  return R::pgamma(1.0/x, alpha, 1.0/beta, false, false);
}


// [[Rcpp::export]]
NumericVector cpp_dinvgamma(
    const NumericVector& x,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), alpha.length(), beta.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    alpha.length(),
    beta.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_invgamma(GETV(x, i), GETV(alpha, i),
                           GETV(beta, i), throw_warning);

  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pinvgamma(
    const NumericVector& x,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), alpha.length(), beta.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    alpha.length(),
    beta.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_invgamma(GETV(x, i), GETV(alpha, i),
                        GETV(beta, i), throw_warning);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}

