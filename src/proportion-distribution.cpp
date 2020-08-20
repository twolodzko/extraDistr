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
*  Re-parametrized beta distribution
*
*  Values:
*  x
*
*  Parameters:
*  0 <= mean <= 1
*  size > 0
*  prior >= 0
*
*/

inline double pdf_prop(double x, double size, double mean, double prior,
                       bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(size) || ISNAN(mean) || ISNAN(prior))
    return x+size+mean+prior;
#endif
  if (size <= 0.0 || mean <= 0.0 || mean >= 1.0 || prior < 0) {
    throw_warning = true;
    return NAN;
  }
  return R::dbeta(x, size*mean+prior, size*(1.0-mean)+prior, false);
}

inline double cdf_prop(double x, double size, double mean, double prior,
                       bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(size) || ISNAN(mean) || ISNAN(prior))
    return x+size+mean+prior;
#endif
  if (size <= 0.0 || mean <= 0.0 || mean >= 1.0 || prior < 0) {
    throw_warning = true;
    return NAN;
  }
  return R::pbeta(x, size*mean+prior, size*(1.0-mean)+prior, true, false);
}

inline double invcdf_prop(double p, double size, double mean, double prior,
                          bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(size) || ISNAN(mean) || ISNAN(prior))
    return p+size+mean+prior;
#endif
  if (size <= 0.0 || mean <= 0.0 || mean >= 1.0 || prior < 0 || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  return R::qbeta(p, size*mean+prior, size*(1.0-mean)+prior, true, false);
}

inline double rng_prop(double size, double mean, double prior,
                       bool& throw_warning) {
  if (ISNAN(size) || ISNAN(mean) || ISNAN(prior) ||
      size <= 0.0 || mean <= 0.0 || mean >= 1.0 || prior < 0) {
    throw_warning = true;
    return NA_REAL;
  }
  return R::rbeta(size*mean+prior, size*(1.0-mean)+prior);
}


// [[Rcpp::export]]
NumericVector cpp_dprop(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& mean,
    const NumericVector& prior,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), size.length(),
                mean.length(), prior.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    size.length(),
    mean.length(),
    prior.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_prop(GETV(x, i), GETV(size, i),
                    GETV(mean, i), GETV(prior, i),
                    throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pprop(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& mean,
    const NumericVector& prior,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), size.length(),
                mean.length(), prior.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    size.length(),
    mean.length(),
    prior.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_prop(GETV(x, i), GETV(size, i),
                    GETV(mean, i), GETV(prior, i),
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
NumericVector cpp_qprop(
    const NumericVector& p,
    const NumericVector& size,
    const NumericVector& mean,
    const NumericVector& prior,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), size.length(),
                mean.length(), prior.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    p.length(),
    size.length(),
    mean.length(),
    prior.length()
  });
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_prop(GETV(pp, i), GETV(size, i),
                       GETV(mean, i), GETV(prior, i),
                       throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rprop(
    const int& n,
    const NumericVector& size,
    const NumericVector& mean,
    const NumericVector& prior
  ) {
  
  if (std::min({size.length(), mean.length(), prior.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_prop(GETV(size, i), GETV(mean, i),
                    GETV(prior, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

