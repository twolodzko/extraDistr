#include <Rcpp.h>

// inline double round0(double x) {
//   return (x < 0.0) ? std::ceil(x) : std::floor(x);
// }

bool isInteger(double x, bool warn) {
  if (ISNAN(x))
    return false;
  if (((x < 0.0) ? std::ceil(x) : std::floor(x)) != x) {
    if (warn) {
      char msg[55];
      std::snprintf(msg, sizeof(msg), "non-integer x = %f", x);
      Rcpp::warning(msg);
    }
    return false;
  }
  return true;
}

// bool anyFinite(Rcpp::NumericVector x) {
//   int n = x.length();
//   for (int i = 0; i < n; i++)
//     if (R_FINITE(x[i]))
//       return true;
//   return false;
// }

double finite_max(Rcpp::NumericVector x) {
  double max_x = 0.0;
  int n = x.length();
  int i = 0;
  do {
    if (R_FINITE(x[i])) {
      max_x = x[i];
      break;
    }
    i++;
  } while (i < n);
  while (i < n) {
    if (R_FINITE(x[i]) && x[i] > max_x) {
      max_x = x[i];
    }
    i++;
  }
  return max_x;
}

// bool allNA(Rcpp::NumericVector x) {
//   int n = x.length();
//   for (int i = 0; i < n; i++)
//     if (!ISNAN(x[i]))
//       return false;
//   return true;
// }

// Random generation

double rng_unif() {
  double u;
  // same as in base R
  do {
    u = R::unif_rand();
  } while (u <= 0.0 || u >= 1.0);
  return u;
}

double rng_sign() {
  double u = rng_unif();
  return (u > 0.5) ? 1.0 : -1.0;
}

