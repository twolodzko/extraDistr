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


inline double logpdf_tpois(double x, double lambda, double a,
                           double b, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(a) || ISNAN(b))
    return x+lambda+a+b;
#endif
  if (lambda < 0.0 || b < a) {
    throw_warning = true;
    return NAN;
  }
  
  if (!isInteger(x) || x < 0.0 || x <= a || x > b || !R_FINITE(x))
    return R_NegInf;
  
  // if (a == 0.0 && b == R_PosInf)
  //   return pow(lambda, x) / (factorial(x) * (exp(lambda) - 1.0));
  
  double pa, pb;
  pa = R::ppois(a, lambda, true, false);
  pb = R::ppois(b, lambda, true, false);
  
  return R::dpois(x, lambda, true) - log(pb-pa);
}

inline double cdf_tpois(double x, double lambda, double a,
                        double b, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(lambda) || ISNAN(a) || ISNAN(b))
    return x+lambda+a+b;
#endif
  if (lambda <= 0.0 || b < a) {
    throw_warning = true;
    return NAN;
  }
  
  if (x < 0.0 || x <= a)
    return 0.0;
  if (x > b || !R_FINITE(x))
    return 1.0;
  
  // if (a == 0.0 && b == R_PosInf)
  //   return R::ppois(x, lambda, true, false) / (1.0 - exp(-lambda));
  
  double pa, pb;
  pa = R::ppois(a, lambda, true, false);
  pb = R::ppois(b, lambda, true, false);

  return (R::ppois(x, lambda, true, false) - pa) / (pb-pa);
}

inline double invcdf_tpois(double p, double lambda, double a,
                           double b, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(lambda) || ISNAN(a) || ISNAN(b))
    return p+lambda+a+b;
#endif
  if (lambda < 0.0 || b < a || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }

  if (p == 0.0)
    return std::max(a, 0.0);
  if (p == 1.0)
    return b;
  
  double pa, pb;
  pa = R::ppois(a, lambda, true, false);
  pb = R::ppois(b, lambda, true, false);
  
  return R::qpois(pa + p*(pb-pa), lambda, true, false);
}

inline double rng_tpois(double lambda, double a, double b,
                        bool& throw_warning) {
  if (ISNAN(lambda) || ISNAN(a) || ISNAN(b) ||
      lambda < 0.0 || b < a) {
    throw_warning = true;
    return NA_REAL;
  }

  double u, pa, pb;
  pa = R::ppois(a, lambda, true, false);
  pb = R::ppois(b, lambda, true, false);
  
  u = R::runif(pa, pb);
  return R::qpois(u, lambda, true, false);
}


// [[Rcpp::export]]
NumericVector cpp_dtpois(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& lower,
    const NumericVector& upper,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), lambda.length(),
                lower.length(), upper.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    lambda.length(),
    lower.length(),
    upper.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_tpois(GETV(x, i), GETV(lambda, i),
                        GETV(lower, i), GETV(upper, i),
                        throw_warning);
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_ptpois(
    const NumericVector& x,
    const NumericVector& lambda,
    const NumericVector& lower,
    const NumericVector& upper,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), lambda.length(),
                lower.length(), upper.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    lambda.length(),
    lower.length(),
    upper.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_tpois(GETV(x, i), GETV(lambda, i),
                     GETV(lower, i), GETV(upper, i),
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
NumericVector cpp_qtpois(
    const NumericVector& p,
    const NumericVector& lambda,
    const NumericVector& lower,
    const NumericVector& upper,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), lambda.length(),
                lower.length(), upper.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    p.length(),
    lambda.length(),
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
    x[i] = invcdf_tpois(GETV(pp, i), GETV(lambda, i),
                        GETV(lower, i), GETV(upper, i),
                        throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rtpois(
    const int& n,
    const NumericVector& lambda,
    const NumericVector& lower,
    const NumericVector& upper
  ) {
  
  if (std::min({lambda.length(), lower.length(), upper.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_tpois(GETV(lambda, i), GETV(lower, i),
                     GETV(upper, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

