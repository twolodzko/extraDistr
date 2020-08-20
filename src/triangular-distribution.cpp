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
*  Triangular distribution
*
*  Values:
*  x
*
*  Parameters:
*  a
*  b > a
*  a <= c <= b
*
*  f(x)    = { (2*(x-a)) / ((b-a)*(c-a))  x < c
*            { 2/(b-a)                    x = c
*            { (2*(b-x)) / ((b-a)*(b-c))  x > c
*  F(x)    = { (x-a)^2 / ((b-a)*(c-a))
*            { 1 - ((b-x)^2 / ((b-a)*(b-c)))
*  F^-1(p) = { a + sqrt(p*(b-a)*(c-a))    p < (c-a)/(b-a)
*            { b - sqrt((1-p)*(b-a)*(b-c));
*/

inline double logpdf_triangular(double x, double a, double b,
                                double c, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(a) || ISNAN(b) || ISNAN(c))
    return x+a+b+c;
#endif
  if (a > c || c > b || a == b) {
    throw_warning = true;
    return NAN;
  }
  if (x < a || x > b) {
    return R_NegInf;
  } else if (x < c) {
    // 2.0*(x-a) / ((b-a)*(c-a));
    return LOG_2F + log(x-a) - log(b-a) - log(c-a);
  } else if (x > c) {
    // 2.0*(b-x) / ((b-a)*(b-c));
    return LOG_2F + log(b-x) - log(b-a) - log(b-c);
  } else {
    // 2.0/(b-a);
    return LOG_2F - log(b-a);
  }
}

inline double cdf_triangular(double x, double a, double b,
                             double c, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(a) || ISNAN(b) || ISNAN(c))
    return x+a+b+c;
#endif
  if (a > c || c > b || a == b) {
    throw_warning = true;
    return NAN;
  }
  if (x < a) {
    return 0.0;
  } else if (x >= b) {
    return 1.0;
  } else if (x <= c) {
    // ((x-a)*(x-a)) / ((b-a)*(c-a));
    return exp( log(x-a) * 2.0 - log(b-a) - log(c-a) );
  } else {
    // 1.0 - (((b-x)*(b-x)) / ((b-a)*(b-c)));
    return 1.0 - exp( log(b-x) * 2.0 - log(b-a) - log(b-c) );
  }
}

inline double invcdf_triangular(double p, double a, double b,
                                double c, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(a) || ISNAN(b) || ISNAN(c))
    return p+a+b+c;
#endif
  if (a > c || c > b || a == b || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  double fc = (c-a)/(b-a);
  if (p < fc)
    return a + sqrt(p*(b-a)*(c-a));
  return b - sqrt((1.0-p)*(b-a)*(b-c));
}

inline double rng_triangular(double a, double b, double c,
                             bool& throw_warning) {
  if (ISNAN(a) || ISNAN(b) || ISNAN(c) ||
      a > c || c > b || a == b) {
    throw_warning = true;
    return NA_REAL;
  }
  double u, v, r, cc;
  r = b - a;
  cc = (c-a)/r;
  u = rng_unif();
  v = rng_unif();
  return ((1.0-cc) * std::min(u, v) + cc * std::max(u, v)) * r + a;
}


// [[Rcpp::export]]
NumericVector cpp_dtriang(
    const NumericVector& x,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), a.length(),
                b.length(), c.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    a.length(),
    b.length(),
    c.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_triangular(GETV(x, i), GETV(a, i),
                          GETV(b, i), GETV(c, i),
                          throw_warning);

  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_ptriang(
    const NumericVector& x,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), a.length(),
                b.length(), c.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    a.length(),
    b.length(),
    c.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_triangular(GETV(x, i), GETV(a, i),
                          GETV(b, i), GETV(c, i),
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
NumericVector cpp_qtriang(
    const NumericVector& p,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), a.length(),
                b.length(), c.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    p.length(),
    a.length(),
    b.length(),
    c.length()
  });
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;

  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;

  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_triangular(GETV(pp, i), GETV(a, i),
                             GETV(b, i), GETV(c, i),
                             throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rtriang(
    const int& n,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c
  ) {
  
  if (std::min({a.length(), b.length(), c.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }

  NumericVector x(n);
  
  bool throw_warning = false;

  for (int i = 0; i < n; i++)
    x[i] = rng_triangular(GETV(a, i), GETV(b, i),
                          GETV(c, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

