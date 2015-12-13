#include <Rcpp.h>
using namespace Rcpp;


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
  return 2 * Phi(x * std::sqrt(2)) - 1;
}

double erfc(double x) {
  return 2 * R::pnorm(x * std::sqrt(2), 0, 1, false, false);
}

double inv_erf(double x) {
  return InvPhi((x+1)/2) / std::sqrt(2);
}

// Factorial

double lfactorial(double x) {
  return R::lgammafn(x + 1);
}

double factorial(double x) {
  return R::gammafn(x + 1);
}

// Random generation for Bernoulli

int rng_bernoulli(double p = 0.5) {
  if (p < 0 || p > 1)
    return NAN;
  double u = R::runif(0, 1);
  if (u <= 1-p)
    return 0;
  else
    return 1;
}
 
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