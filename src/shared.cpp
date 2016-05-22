#include <Rcpp.h>
#include "const.h"
using namespace Rcpp;

// Basic functions

bool tol_equal(double x, double y) {
  return std::abs(x - y) <= MIN_DIFF_EPS;
}

bool isInteger(double x) {
  if (floor(x) != x) {
    char msg[55];
    std::snprintf(msg, sizeof(msg), "non-integer x = %f", x);
    Rcpp::warning(msg);
    return false;
  }
  return true;
}

// Dealing with Inf

bool anyFinite(NumericVector x) {
  int n = x.length();
  for (int i = 0; i < n; i++)
    if (!std::isinf(x[i]))
      return true;
  return false;
}

double finite_max(NumericVector x) {
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
  return R::dnorm(x, 0, 1, false);
}

double Phi(double x) {
  return R::pnorm(x, 0, 1, true, false);
}

double InvPhi(double x) {
  return R::qnorm(x, 0, 1, true, false);
}

// Error function

double erf(double x) {
  return 2 * Phi(x * sqrt(2)) - 1;
}

double erfc(double x) {
  return 2 * R::pnorm(x * sqrt(2), 0, 1, false, false);
}

double inv_erf(double x) {
  return InvPhi((x+1)/2) / sqrt(2);
}

// Factorial

double lfactorial(double x) {
  return R::lgammafn(x + 1);
}

double factorial(double x) {
  return R::gammafn(x + 1);
}

// Random generation for Bernoulli

double rng_bernoulli(double p = 0.5) {
  if (p < 0 || p > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double u = R::runif(0, 1);
  if (u <= 1-p)
    return 0;
  else
    return 1;
}
 

//  // Multivariate gamma function
//  
// double lmvgammafn(double p, double a) {
//   double prod_gamma = 0;
//   for (int j = 1; j <= p; j++)
//     prod_gamma += R::lgammafn(a+(1-j)/2);
//   return M_PI * (p*(p-1)/4) + prod_gamma;
// }
//  
 
// /*
//  * Incomplete gamma function
//  * 
//  * pgamma(x, a, ..) currently requires a > 0, whereas the
//  * incomplete gamma function is also defined for negative a
//  * 
//  */
// 
// double inc_gamma(double x, double a) {
//   return R::pgamma(x, a, 1, true, false) * R::gammafn(a);
// }
// 
// double inc_ugamma(double x, double a) {
//   return R::pgamma(x, a, 1, false, false) * R::gammafn(a);
// }
// 
// double log_inc_gamma(double x, double a) {
//   return R::pgamma(x, a, 1, true, true) + R::lgammafn(a);
// }
// 
// double log_inc_ugamma(double x, double a) {
//   return R::pgamma(x, a, 1, false, true) + R::lgammafn(a);
// }