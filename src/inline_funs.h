
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

inline double to_dbl(long int x) {
  return static_cast<double>(x);
}

inline long int to_int(double x) {
  if (R_FINITE(x) && x > std::numeric_limits<long int>::max())
    Rcpp::stop("reached largest integer which can be represented as <long int>");
  return static_cast<long int>(x);
}


#endif

