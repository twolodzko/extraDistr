#include <Rcpp.h>
#include "shared.h"
// [[Rcpp::plugins(cpp11)]]

using Rcpp::NumericVector;


// [[Rcpp::export]]
NumericVector cpp_rsign(
    const int& n
  ) {
  
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_sign();
  
  return x;
}

