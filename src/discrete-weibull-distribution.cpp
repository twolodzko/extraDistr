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
*  Discrete Weibull distribution
*
*  Values:
*  x >= 0
*
*  Parameters:
*  0 < q < 1
*  beta
*
*  f(x)    = q^x^beta - q^(x+1)^beta
*  F(x)    = 1-q^(x+1)^beta
*  F^-1(p) = ceiling(pow(log(1-p)/log(q), 1/beta) - 1)
*
*  Nakagawa and Osaki (1975), "The Discrete Weibull Distribution",
*  IEEE Transactions on Reliability, R-24, pp. 300-301.
*
*/

inline double pdf_dweibull(double x, double q, double beta,
                           bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(q) || ISNAN(beta))
    return x+q+beta;
#endif
  if (q <= 0.0 || q >= 1.0 || beta <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (!isInteger(x) || x < 0.0)
    return 0.0;
  return pow(q, pow(x, beta)) - pow(q, pow(x+1.0, beta));
}

inline double cdf_dweibull(double x, double q, double beta,
                           bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(q) || ISNAN(beta))
    return x+q+beta;
#endif
  if (q <= 0.0 || q >= 1.0 || beta <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  // 1.0 - pow(q, pow(x+1.0, beta))
  return 1.0 - exp(log(q) * exp(log1p(x) * beta));
}

inline double invcdf_dweibull(double p, double q, double beta,
                              bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(q) || ISNAN(beta))
    return p+q+beta;
#endif
  if (q <= 0.0 || q >= 1.0 || beta <= 0.0 || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  if (p == 0.0)
    return 0.0;
  return ceil(pow(log(1.0 - p)/log(q), 1.0/beta) - 1.0);
}

inline double rng_dweibull(double q, double beta,
                           bool& throw_warning) {
  if (ISNAN(q) || ISNAN(beta) || q <= 0.0 || q >= 1.0 ||
      beta <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  double u = rng_unif();
  return ceil(pow(log(u)/log(q), 1.0/beta) - 1.0);
}


// [[Rcpp::export]]
NumericVector cpp_ddweibull(
    const NumericVector& x,
    const NumericVector& q,
    const NumericVector& beta,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), q.length(), beta.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    q.length(),
    beta.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_dweibull(GETV(x, i), GETV(q, i),
                        GETV(beta, i), throw_warning);

  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pdweibull(
    const NumericVector& x,
    const NumericVector& q,
    const NumericVector& beta,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), q.length(), beta.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    q.length(),
    beta.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_dweibull(GETV(x, i), GETV(q, i),
                        GETV(beta, i), throw_warning);

  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qdweibull(
    const NumericVector& p,
    const NumericVector& q,
    const NumericVector& beta,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), q.length(), beta.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    p.length(),
    q.length(),
    beta.length()
  });
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;

  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;

  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_dweibull(GETV(pp, i), GETV(q, i),
                           GETV(beta, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rdweibull(
    const int& n,
    const NumericVector& q,
    const NumericVector& beta
  ) {
  
  if (std::min({q.length(), beta.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }

  NumericVector x(n);
  
  bool throw_warning = false;

  for (int i = 0; i < n; i++)
    x[i] = rng_dweibull(GETV(q, i), GETV(beta, i),
                        throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

