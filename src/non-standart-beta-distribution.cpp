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
*  Non-standard beta distribution
*
*  Values:
*  x
*
*  Parameters:
*  0 <= beta <= 1
*  alpha > 0
*  lower < upper
*
*/

double pdf_nsbeta(double x, double alpha, double beta, double l,
                  double u, bool log_p, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(alpha) || ISNAN(beta) || ISNAN(l) || ISNAN(u))
    return x+alpha+beta+l+u;
#endif
  if (l >= u || alpha < 0.0 || beta < 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double r = u-l;
  double p = R::dbeta((x-l)/r, alpha, beta, log_p);
  if (log_p) 
    return p-log(r);
  else
    return p/r;
}

double cdf_nsbeta(double x, double alpha, double beta, double l,
                  double u, bool lower_tail, bool log_p, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(alpha) || ISNAN(beta) || ISNAN(l) || ISNAN(u))
    return x+alpha+beta+l+u;
#endif
  if (l >= u || alpha < 0.0 || beta < 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return R::pbeta((x-l)/(u-l), alpha, beta, lower_tail, log_p);
}

double invcdf_nsbeta(double p, double alpha, double beta, double l,
                     double u, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(alpha) || ISNAN(beta) || ISNAN(l) || ISNAN(u))
    return p+alpha+beta+l+u;
#endif
  if (l >= u || alpha < 0.0 || beta < 0.0 || !VALID_PROB(p)) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return R::qbeta(p, alpha, beta, true, false) * (u-l) + l;
}

double rng_nsbeta(double alpha, double beta, double l, double u,
                  bool& throw_warning) {
  if (ISNAN(alpha) || ISNAN(beta) || ISNAN(l) || ISNAN(u) ||
      l >= u || alpha < 0.0 || beta < 0.0) {
    Rcpp::warning("NAs produced");
    return NA_REAL;
  }
  return R::rbeta(alpha, beta) * (u-l) + l;
}


// [[Rcpp::export]]
NumericVector cpp_dnsbeta(
    const NumericVector& x,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& lower,
    const NumericVector& upper,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), alpha.length(),
                beta.length(), lower.length(),
                upper.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    alpha.length(),
    beta.length(),
    lower.length(),
    upper.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_nsbeta(GETV(x, i), GETV(alpha, i),
                      GETV(beta, i), GETV(lower, i),
                      GETV(upper, i), log_prob, throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pnsbeta(
    const NumericVector& x,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& lower,
    const NumericVector& upper,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), alpha.length(),
                beta.length(), lower.length(),
                upper.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    alpha.length(),
    beta.length(),
    lower.length(),
    upper.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_nsbeta(GETV(x, i), GETV(alpha, i),
                      GETV(beta, i), GETV(lower, i),
                      GETV(upper, i), lower_tail,
                      log_prob, throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qnsbeta(
    const NumericVector& p,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& lower,
    const NumericVector& upper,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), alpha.length(),
                beta.length(), lower.length(),
                upper.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    p.length(),
    alpha.length(),
    beta.length(),
    lower.length(),
    upper.length()
  });
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_nsbeta(GETV(pp, i), GETV(alpha, i),
                         GETV(beta, i), GETV(lower, i),
                         GETV(upper, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rnsbeta(
    const int& n,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& lower,
    const NumericVector& upper
  ) {
  
  if (std::min({alpha.length(), beta.length(),
                lower.length(), upper.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_nsbeta(GETV(alpha, i), GETV(beta, i),
                      GETV(lower, i), GETV(upper, i),
                      throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

