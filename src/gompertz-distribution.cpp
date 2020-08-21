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
*  Gompertz distribution
*
*  Values:
*  x >= 0
*
*  Parameters:
*  a > 0
*  b > 0
*
*  f(x)    = a*exp(b*x - a/b * (exp(bx)-1))
*  F(x)    = 1-exp(-a/b * (exp(bx)-1))
*  F^-1(p) = 1/b * log(1 - b/a * log(1-p))
*
* References:
*
* Lenart, A. (2012). The Gompertz distribution and Maximum Likelihood Estimation
* of its parameters - a revision. MPIDR WORKING PAPER WP 2012-008.
*
*/


inline double logpdf_gompertz(double x, double a, double b,
                              bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(a) || ISNAN(b))
    return x+a+b;
#endif
  if (a <= 0.0 || b <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0 || !R_FINITE(x))
    return R_NegInf;
  // a * exp(b*x - a/b * (exp(b*x) - 1.0));
  return log(a) + (b*x - a/b * (exp(b*x) - 1.0));
}

inline double cdf_gompertz(double x, double a, double b,
                           bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(a) || ISNAN(b))
    return x+a+b;
#endif
  if (a <= 0.0 || b <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  if (!R_FINITE(x))
    return 1.0;
  return 1.0 - exp(-a/b * (exp(b*x) - 1.0));
}

inline double invcdf_gompertz(double p, double a, double b,
                              bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(a) || ISNAN(b))
    return p+a+b;
#endif
  if (a <= 0.0 || b <= 0.0 || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  return log(1.0 - b/a * log(1.0-p)) / b;
}

inline double rng_gompertz(double a, double b, bool& throw_warning) {
  if (ISNAN(a) || ISNAN(b) || a <= 0.0 || b <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  double u = rng_unif();
  return log(1.0 - b/a * log(u)) / b;
}


// [[Rcpp::export]]
NumericVector cpp_dgompertz(
    const NumericVector& x,
    const NumericVector& a,
    const NumericVector& b,
    bool log_prob = false
  ) {
  
  if (std::min({x.length(), a.length(), b.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    a.length(),
    b.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_gompertz(GETV(x, i), GETV(a, i),
                           GETV(b, i), throw_warning);

  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pgompertz(
    const NumericVector& x,
    const NumericVector& a,
    const NumericVector& b,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), a.length(), b.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    a.length(),
    b.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_gompertz(GETV(x, i), GETV(a, i),
                        GETV(b, i), throw_warning);

  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qgompertz(
    const NumericVector& p,
    const NumericVector& a,
    const NumericVector& b,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), a.length(), b.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    p.length(),
    a.length(),
    b.length()
  });
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_gompertz(GETV(pp, i), GETV(a, i),
                           GETV(b, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rgompertz(
    const int& n,
    const NumericVector& a,
    const NumericVector& b
  ) {
  
  if (std::min({a.length(), b.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }

  NumericVector x(n);
  
  bool throw_warning = false;

  for (int i = 0; i < n; i++)
    x[i] = rng_gompertz(GETV(a, i), GETV(b, i),
                        throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

