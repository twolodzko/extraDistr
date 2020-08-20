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
*  Kumaraswamy distribution
*
*  Values:
*  x in [0, 1]
*
*  Parameters:
*  a > 0
*  b > 0
*
*  f(x)    = a*b*x^{a-1}*(1-x^a)^{b-1}
*  F(x)    = 1-(1-x^a)^b
*  F^-1(p) = 1-(1-p^{1/b})^{1/a}
*
*/

inline double pdf_kumar(double x, double a, double b,
                        bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(a) || ISNAN(b))
    return x+a+b;
#endif
  if (a <= 0.0 || b <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0 || x > 1.0)
    return 0.0;
  // is the support [0,1] or (0,1) ?
  // log(a) + log(b) + log(x)*(a-1.0) + log1p(-pow(x, a))*(b-1.0);
  return a*b * pow(x, a-1.0) * pow(1.0-pow(x, a), b-1.0);
}

inline double cdf_kumar(double x, double a, double b,
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
  if (x >= 1.0)
    return 1.0;
  return 1.0 - pow(1.0 - pow(x, a), b);
}

inline double invcdf_kumar(double p, double a, double b,
                           bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(a) || ISNAN(b))
    return p+a+b;
#endif
  if (a <= 0.0 || b <= 0.0 || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  return pow(1.0 - pow(1.0 - p, 1.0/b), 1.0/a);
}

inline double rng_kumar(double a, double b, bool& throw_warning) {
  if (ISNAN(a) || ISNAN(b) || a <= 0.0 || b <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  double u = rng_unif();
  return pow(1.0 - pow(u, 1.0/b), 1.0/a);
}


// [[Rcpp::export]]
NumericVector cpp_dkumar(
    const NumericVector& x,
    const NumericVector& a,
    const NumericVector& b,
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
    p[i] = pdf_kumar(GETV(x, i), GETV(a, i),
                     GETV(b, i), throw_warning);

  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pkumar(
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
    p[i] = cdf_kumar(GETV(x, i), GETV(a, i),
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
NumericVector cpp_qkumar(
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
    q[i] = invcdf_kumar(GETV(pp, i), GETV(a, i),
                        GETV(b, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rkumar(
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
    x[i] = rng_kumar(GETV(a, i), GETV(b, i),
                     throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

