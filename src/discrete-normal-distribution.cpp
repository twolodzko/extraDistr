#include <Rcpp.h>
#include "shared.h"

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


double pmf_dnorm(double x, double mu,
                 double sigma, bool& throw_warning) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
    return x+mu+sigma;
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
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pmf_dnorm(x[i % dims[0]], mu[i % dims[1]],
                     sigma[i % dims[2]], throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}

