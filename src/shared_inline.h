
#ifndef EDCPP_INLINEFUNS_H
#define EDCPP_INLINEFUNS_H

#include "shared.h"
#include <Rcpp.h>


inline bool tol_equal(double x, double y) {
  return std::abs(x - y) < MIN_DIFF_EPS;
}

inline double phi(double x) {
  return R::dnorm(x, 0.0, 1.0, false);
}

inline double lphi(double x) {
  return R::dnorm(x, 0.0, 1.0, true);
}

inline double Phi(double x) {
  return R::pnorm(x, 0.0, 1.0, true, false);
}

inline double InvPhi(double x) {
  return R::qnorm(x, 0.0, 1.0, true, false);
}

inline double factorial(double x) {
  return R::gammafn(x + 1.0);
}

inline double lfactorial(double x) {
  return R::lgammafn(x + 1.0);
}

inline double rng_sign() {
  double u = rng_unif();
  return (u > 0.5) ? 1.0 : -1.0;
}

inline bool is_large_int(double x) {
  if (x > std::numeric_limits<int>::max())
    return true;
  return false;
}

inline double to_dbl(int x) {
  return static_cast<double>(x);
}

inline int to_pos_int(double x) {
  if (x < 0.0 || ISNAN(x))
    Rcpp::stop("value cannot be coerced to integer");
  if (is_large_int(x))
    Rcpp::stop("value out of integer range");
  return static_cast<int>(x);
}

inline double trunc_p(double x) {
  return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x); 
}


#endif

