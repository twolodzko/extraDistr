#include <Rcpp.h>
#include "shared.h"
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]

using std::log;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;


/*
 * 
 *
 * Random generation from categorical distribution using Gumbel-max trick
 * 
 * Random generation from categorical distribution given log-probabilities.
 * 
 * Given vector of log-probabilities gamma[1], ..., gamma[k], where each
 * gamma[i] = log(prob[i]) + constant, then it is possible to taking 
 * 
 * x = argmax[i]( gamma[i] + g[i] )
 * 
 * where g[1], ..., g[k] is a sample from standard Gumbel distribution.
 * 
 * 
 * References:
 * 
 * Maddison, C. J., Tarlow, D., & Minka, T. (2014). A* sampling.
 * [In:] Advances in Neural Information Processing Systems (pp. 3086-3094).
 * 
 * 
 */


// [[Rcpp::export]]
NumericVector cpp_rcatlp(
    const int& n,
    const NumericMatrix& log_prob
  ) {
  
  if (log_prob.length() < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  int k = log_prob.ncol();
  double u, glp, max_val;
  int jj;
  
  bool throw_warning = false;
  bool wrong_prob = false;
  
  for (int i = 0; i < n; i++) {
    
    max_val = -INFINITY;
    jj = 0;
    
    for (int j = 0; j < k; j++) {
      
      wrong_prob = false;
      
      if (ISNAN(GETM(log_prob, i, j))) {
        throw_warning = wrong_prob = true;
        break;
      }
      
      u = R::exp_rand(); // -log(rng_unif())
      glp = -log(u) + GETM(log_prob, i, j); 
      if (glp > max_val) {
        max_val = glp;
        jj = j;
      }
      
    }
    
    if (wrong_prob) {
      x[i] = NA_REAL;
    } else {
      x[i] = static_cast<double>(jj);
    }
  }
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

