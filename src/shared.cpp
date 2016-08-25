#include <Rcpp.h>
#include "const.h"


// Basic functions

bool tol_equal(double x, double y) {
  return std::abs(x - y) < MIN_DIFF_EPS;
}

bool isInteger(double x) {
  if (std::floor(x) != x) {
    char msg[55];
    std::snprintf(msg, sizeof(msg), "non-integer x = %f", x);
    Rcpp::warning(msg);
    return false;
  }
  return true;
}

// Dealing with Inf

bool anyFinite(Rcpp::NumericVector x) {
  int n = x.length();
  for (int i = 0; i < n; i++)
    if (!std::isinf(x[i]))
      return true;
  return false;
}

double finite_max(Rcpp::NumericVector x) {
  double max_x = -INFINITY;
  int n = x.length();
  for (int i = 0; i < n; i++) {
    if (!std::isinf(x[i]) && x[i] > max_x)
      max_x = x[i];
  }
  return max_x;
}

// Standard normal

double phi(double x) {
  return R::dnorm(x, 0.0, 1.0, false);
}

double Phi(double x) {
  return R::pnorm(x, 0.0, 1.0, true, false);
}

double InvPhi(double x) {
  return R::qnorm(x, 0.0, 1.0, true, false);
}

// Factorial

double factorial(double x) {
  return R::gammafn(x + 1.0);
}

double lfactorial(double x) {
  return R::lgammafn(x + 1.0);
}

// Random generation

double rng_unif() {
  double u;
  // same as in base R
  do {
    u = R::unif_rand();
  } while (u <= 0.0 || u >= 1.0);
  return u;
}

double rng_bern(double p) {
  if (p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double u = rng_unif();
  return (u > p) ? 0.0 : 1.0;
}

double rng_sign() {
  double u = rng_unif();
  return (u > 0.5) ? 1.0 : -1.0;
}

// Checking parameters

Rcpp::NumericMatrix normalize_prob(const Rcpp::NumericMatrix& prob) {
  
  int n = prob.nrow();
  int k = prob.ncol();
  double p_tot;
  bool wrong_param;
  Rcpp::NumericMatrix p = Rcpp::clone(prob);
  
  for (int i = 0; i < n; i++) {
    wrong_param = false;
    p_tot = 0.0;
    
    for (int j = 0; j < k; j++) {
      if (p(i, j) < 0.0 || ISNAN(p(i, j))) {
        wrong_param = true;
        break;
      }
      p_tot += p(i, j);
    }
    
    if (wrong_param) {
      Rcpp::warning("NaNs produced");
      for (int j = 0; j < k; j++)
        p(i, j) = NAN;
    } else if (p_tot > 1.0) {
      for (int j = 0; j < k; j++)
        p(i, j) /= p_tot;
    }
  }
  
  return p;
}

