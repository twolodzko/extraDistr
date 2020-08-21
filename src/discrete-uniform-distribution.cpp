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
 * Discrete uniform distribution
 * 
 * Values:
 * a <= x <= b
 * 
 * f(x) = 1/(b-a+1)
 * F(x) = (floor(x)-a+1)/b-a+1
 *  
 */


inline double pmf_dunif(double x, double min, double max,
                        bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(min) || ISNAN(max))
    return x+min+max;
#endif
  if (min > max || !R_FINITE(min) || !R_FINITE(max) ||
      !isInteger(min, false) || !isInteger(max, false)) {
    throw_warning = true;
    return NAN;
  }
  if (x < min || x > max || !isInteger(x))
    return 0.0;
  return 1.0/(max-min+1.0);
}


inline double cdf_dunif(double x, double min, double max,
                        bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(min) || ISNAN(max))
    return x+min+max;
#endif
  if (min > max || !R_FINITE(min) || !R_FINITE(max) ||
      !isInteger(min, false) || !isInteger(max, false)) {
    throw_warning = true;
    return NAN;
  }
  if (x < min)
    return 0.0;
  else if (x >= max)
    return 1.0;
  return (floor(x)-min+1.0)/(max-min+1.0);
}

inline double invcdf_dunif(double p, double min, double max,
                           bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(min) || ISNAN(max))
    return p+min+max;
#endif
  if (min > max || !R_FINITE(min) || !R_FINITE(max) ||
      !isInteger(min, false) || !isInteger(max, false) ||
      !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  if (p == 0 || min == max)
    return min;
  return ceil( p*(max-min+1.0)+min-1.0 );
}

inline double rng_dunif(double min, double max, bool& throw_warning) {
  if (ISNAN(min) || ISNAN(max) ||
      min > max || !R_FINITE(min) || !R_FINITE(max) ||
      !isInteger(min, false) || !isInteger(max, false)) {
    throw_warning = true;
    return NA_REAL;
  }
  if (min == max)
    return min;
  return ceil(R::runif(min - 1.0, max));
}


// [[Rcpp::export]]
NumericVector cpp_ddunif(
    const NumericVector& x,
    const NumericVector& min,
    const NumericVector& max,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), min.length(), max.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    min.length(),
    max.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pmf_dunif(GETV(x, i), GETV(min, i),
                     GETV(max, i), throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pdunif(
    const NumericVector& x,
    const NumericVector& min,
    const NumericVector& max,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), min.length(), max.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    min.length(),
    max.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_dunif(GETV(x, i), GETV(min, i),
                     GETV(max, i), throw_warning);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qdunif(
    const NumericVector& p,
    const NumericVector& min,
    const NumericVector& max,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), min.length(), max.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    p.length(),
    min.length(),
    max.length()
  });
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_dunif(GETV(pp, i), GETV(min, i),
                        GETV(max, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rdunif(
    const int& n,
    const NumericVector& min,
    const NumericVector& max
  ) {
  
  if (std::min({min.length(), max.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_dunif(GETV(min, i), GETV(max, i),
                     throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

