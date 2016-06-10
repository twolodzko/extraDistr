#include <Rcpp.h>
#include "shared.h"

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using std::sin;
using std::cos;
using std::tan;
using std::atan;
using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;


/*
 * Skellam distribution
 * 
 * mu1 >= 0
 * mu2 >= 0
 * 
 */

double pmf_skellam(double x, double mu1, double mu2) {
  if (mu1 < 0.0 || mu2 <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(x))
    return 0.0;
  return exp(-(mu1+mu2)) * pow(mu1/mu2, x/2.0) * R::bessel_i(2*sqrt(mu1*mu2), x, 1.0);
}

double rng_skellam(double mu1, double mu2) {
  if (mu1 < 0.0 || mu2 <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return R::rpois(mu1) - R::rpois(mu2);
}



// [[Rcpp::export]]
NumericVector cpp_dskellam(
    const NumericVector& x,
    const NumericVector& mu1,
    const NumericVector& mu2,
    bool log_prob = false
  ) {
  
  int nx = x.length();
  int na = mu1.length();
  int nb = mu2.length();
  int Nmax = Rcpp::max(IntegerVector::create(nx, na, nb));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pmf_skellam(x[i % nx], mu1[i % na], mu2[i % nb]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rskellam(
    const int n,
    const NumericVector& mu1,
    const NumericVector& mu2
  ) {
  
  int na = mu1.length();
  int nb = mu2.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_skellam(mu1[i % na], mu2[i % nb]);
  
  return x;
}

