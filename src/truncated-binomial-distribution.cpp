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


inline double logpdf_tbinom(double x, double size, double prob, double a,
                            double b, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(size) || ISNAN(prob) || ISNAN(a) || ISNAN(b))
    return x+size+prob+a+b;
#endif
  if (size < 0.0 || !VALID_PROB(prob) || b < a || !isInteger(size, false)) {
    throw_warning = true;
    return NAN;
  }
  
  if (!isInteger(x) || x < 0.0 || x <= a || x > b || x > size)
    return R_NegInf;
  
  double pa, pb;
  pa = R::pbinom(a, size, prob, true, false);
  pb = R::pbinom(b, size, prob, true, false);
  
  return R::dbinom(x, size, prob, true) - log(pb-pa);
}

inline double cdf_tbinom(double x, double size, double prob, double a,
                         double b, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(size) || ISNAN(prob) || ISNAN(a) || ISNAN(b))
    return x+size+prob+a+b;
#endif
  if (size < 0.0 || !VALID_PROB(prob) || b < a ||
      !isInteger(size, false)) {
    throw_warning = true;
    return NAN;
  }
  
  if (x < 0.0 || x <= a)
    return 0.0;
  if (x > b || x >= size)
    return 1.0;
  
  double pa, pb;
  pa = R::pbinom(a, size, prob, true, false);
  pb = R::pbinom(b, size, prob, true, false);
  
  return (R::pbinom(x, size, prob, true, false) - pa) / (pb-pa);
}

inline double invcdf_tbinom(double p, double size, double prob,
                            double a, double b, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(size) || ISNAN(prob) || ISNAN(a) || ISNAN(b))
    return p+size+prob+a+b;
#endif
  if (size < 0.0 || !VALID_PROB(prob) || b < a ||
      !isInteger(size, false) || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  
  if (p == 0.0)
    return std::max(a, 0.0);
  if (p == 1.0)
    return std::min(size, b);
  
  double pa, pb;
  pa = R::pbinom(a, size, prob, true, false);
  pb = R::pbinom(b, size, prob, true, false);
  
  return R::qbinom(pa + p*(pb-pa), size, prob, true, false);
}

inline double rng_tbinom(double size, double prob, double a,
                         double b, bool& throw_warning) {
  if (ISNAN(size) || ISNAN(prob) || ISNAN(a) || ISNAN(b) ||
      size < 0.0 || !VALID_PROB(prob) || b < a ||
      !isInteger(size, false)) {
    throw_warning = true;
    return NA_REAL;
  }
  
  double u, pa, pb;
  pa = R::pbinom(a, size, prob, true, false);
  pb = R::pbinom(b, size, prob, true, false);
  
  u = R::runif(pa, pb);
  return R::qbinom(u, size, prob, true, false);
}


// [[Rcpp::export]]
NumericVector cpp_dtbinom(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& lower,
    const NumericVector& upper,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), size.length(), prob.length(),
                lower.length(), upper.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    size.length(),
    prob.length(),
    lower.length(),
    upper.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_tbinom(GETV(x, i), GETV(size, i),
                         GETV(prob, i), GETV(lower, i),
                         GETV(upper, i), throw_warning);
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_ptbinom(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& lower,
    const NumericVector& upper,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), size.length(), prob.length(),
                lower.length(), upper.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    size.length(),
    prob.length(),
    lower.length(),
    upper.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_tbinom(GETV(x, i), GETV(size, i),
                      GETV(prob, i), GETV(lower, i),
                      GETV(upper, i), throw_warning);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qtbinom(
    const NumericVector& p,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& lower,
    const NumericVector& upper,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), size.length(), prob.length(),
                lower.length(), upper.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    p.length(),
    size.length(),
    prob.length(),
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
    x[i] = invcdf_tbinom(GETV(pp, i), GETV(size, i),
                         GETV(prob, i), GETV(lower, i),
                         GETV(upper, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rtbinom(
    const int& n,
    const NumericVector& size,
    const NumericVector& prob,
    const NumericVector& lower,
    const NumericVector& upper
  ) {
  
  if (std::min({size.length(), prob.length(),
                lower.length(), upper.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_tbinom(GETV(size, i), GETV(prob, i),
                      GETV(lower, i), GETV(upper, i),
                      throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

