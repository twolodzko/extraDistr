#include <Rcpp.h>
#include "const.h"
#include "shared.h"

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;


/*
* Multinomial distribution
* 
* x[i]       number of values of i-th category drawn
* n = sum(x) total number of draws
* p[i]       probability of drawing i-th category value
*  
* f(x) = n!/prod(x[i]!) * prod(p[i]^x[i])
*
*/


// [[Rcpp::export]]
NumericVector cpp_dmnom(
    const NumericMatrix& x,
    const NumericVector& size,
    const NumericMatrix& prob,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.nrow());
  dims.push_back(size.length());
  dims.push_back(prob.nrow());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  int m = x.ncol();
  int k = prob.ncol();
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  if (m != k)
    Rcpp::stop("number of columns in x does not equal number of columns in prob");
  
  double n_fac, prod_xfac, prod_pow_px, sum_x, p_tot;
  bool wrong_param, wrong_x;
  
  for (int i = 0; i < Nmax; i++) {
    
    sum_x = 0.0;
    p_tot = 0.0;
    wrong_param = false;
    wrong_x = false;
    
    for (int j = 0; j < k; j++) {
      if (prob(i % dims[2], j) < 0.0)
        wrong_param = true;
      p_tot += prob(i % dims[2], j);
      sum_x += x(i % dims[0], j);
    }
    
    if (ISNAN(p_tot + sum_x + size[i % dims[1]])) {
      p[i] = p_tot + sum_x + size[i % dims[1]];
      continue;
    } 
    
    if (wrong_param || size[i % dims[1]] < 0.0 ||
        !isInteger(size[i % dims[1]], false)) {
      throw_warning = true;
      p[i] = NAN; 
      continue;
    }
    
    n_fac = lfactorial(size[i % dims[1]]);
    prod_xfac = 0.0;
    prod_pow_px = 0.0;
    
    for (int j = 0; j < k; j++) {
      if (x(i % dims[0], j) < 0.0 || !isInteger(x(i % dims[0], j))) {
        wrong_x = true;
      } else {
        prod_xfac += lfactorial(x(i % dims[0], j));
        prod_pow_px += log(prob(i % dims[2], j) / p_tot) * x(i % dims[0], j);
      }
    }
    
    if (wrong_x || sum_x < 0.0 || sum_x != size[i % dims[1]]) {
      p[i] = R_NegInf;
    } else {
      p[i] = n_fac - prod_xfac + prod_pow_px;
    }
  }
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericMatrix cpp_rmnom(
    const int& n,
    const NumericVector& size,
    const NumericMatrix& prob
  ) {
  
  std::vector<int> dims;
  dims.push_back(size.length());
  dims.push_back(prob.nrow());
  int k = prob.ncol();
  bool wrong_values;
  double p_tmp, size_left, sum_p, p_tot;
  
  NumericMatrix x(n, k);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++) {
    
    size_left = size[i % dims[0]];
    sum_p = 1.0;
    p_tot = 0.0;
    wrong_values = false;
    
    // TODO:
    // sort prob(i,_) first?
    
    for (int j = 0; j < k; j++) {
      if (prob(i % dims[1], j) < 0.0) {
        wrong_values = true;
        break;
      }
      p_tot += prob(i % dims[1], j);
    }
    
    if (wrong_values || ISNAN(p_tot + size[i % dims[0]]) ||
        size[i % dims[0]] < 0.0 || !isInteger(size[i % dims[0]], false)) {
      throw_warning = true;
      for (int j = 0; j < k; j++)
        x(i, j) = NA_REAL;
      continue;
    }

    for (int j = 0; j < k-1; j++) {
      p_tmp = prob(i % dims[1], j)/p_tot;
      x(i, j) = R::rbinom(size_left, p_tmp/sum_p);
      size_left -= x(i, j);
      sum_p -= p_tmp;
    }
    
    x(i, k-1) = size_left;
    
  }
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

