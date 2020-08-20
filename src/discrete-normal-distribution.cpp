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


inline double pmf_dnorm(double x, double mu, double sigma,
                        bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
    return x+mu+sigma;
#endif
  if (sigma <= 0.0) {
    throw_warning = true;
    return NAN;
  }
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
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), mu.length(), sigma.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    mu.length(),
    sigma.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pmf_dnorm(GETV(x, i), GETV(mu, i),
                     GETV(sigma, i), throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}

