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
* Discrete normal distribution
* 
* Values:
* x
* 
* Parameters
* mu
* sigma > 0
*  
*/


double pmf_dnorm(double x, double mu, double sigma) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
    return NAN;
  // if (sigma <= 0.0) {
  //   Rcpp::warning("NaNs produced");
  //   return NAN;
  // }
  if (!isInteger(x))
    return 0.0;
  return R::pnorm(x+1.0, mu, sigma, true, false) -
         R::pnorm(x, mu, sigma, true, false);
}


// [[Rcpp::export]]
NumericVector cpp_ddnorm(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    bool log_prob = false
  ) {
  
  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns));
  NumericVector p(Nmax);
  NumericVector sigma_n = positive_or_nan(sigma);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pmf_dnorm(x[i % n], mu[i % nm], sigma_n[i % ns]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}

