#include <Rcpp.h>
#include "shared.h"

bool isInteger(double x, bool warn) {
  if (ISNAN(x))
    return false;
  if (((x < 0.0) ? std::ceil(x) : std::floor(x)) != x) {
    if (warn) {
      char msg[55];
      std::snprintf(msg, sizeof(msg), "non-integer: %f", x);
      Rcpp::warning(msg);
    }
    return false;
  }
  return true;
}

double finite_max_int(const Rcpp::NumericVector& x) {
  double max_x = 0.0;
  int n = x.length();
  int i = 0;
  do {
    if (x[i] > 0.0 && !is_large_int(x[i])) {
      max_x = x[i];
      break;
    }
    i++;
  } while (i < n);
  while (i < n) {
    if (x[i] > max_x && !is_large_int(x[i])) {
      max_x = x[i];
    }
    i++;
  }
  return max_x;
}

double rng_unif() {
  double u;
  // same as in base R
  do {
    u = R::unif_rand();
  } while (u <= 0.0 || u >= 1.0);
  return u;
}

