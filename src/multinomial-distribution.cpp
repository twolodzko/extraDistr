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
using std::sin;
using std::cos;
using std::tan;
using std::atan;
using Rcpp::IntegerVector;
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
    bool log_prob = false
  ) {
  
  int n = x.nrow();
  int m = x.ncol();
  int k = prob.ncol();
  int np = prob.nrow();
  int ns = size.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, ns, np));
  NumericVector p(Nmax);
  
  if (m != k)
    Rcpp::stop("Number of columns in 'x' does not equal number of columns in 'prob'.");
  
  double n_fac, prod_xfac, prod_pow_px, sum_x, p_tot;
  bool wrong_param, wrong_x, missings;
  
  for (int i = 0; i < Nmax; i++) {
    
    n_fac = lfactorial(size[i % ns]);
    prod_xfac = 0.0;
    prod_pow_px = 0.0;
    
    sum_x = 0.0;
    p_tot = 0.0;
    wrong_param = false;
    wrong_x = false;
    missings = false;
    
    if (ISNAN(size[i % ns]))
      missings = true;
    
    for (int j = 0; j < k; j++) {
      if (ISNAN(prob(i % np, j)) || ISNAN(x(i % n, j))) {
        missings = true;
        break;
      }
    }
    
    if (missings) {
      p[i] = NA_REAL;
      continue;
    } 
    
    if (size[i % ns] < 0.0 || floor(size[i % ns]) != size[i % ns]) {
      wrong_param = true;
    } else {
      
      for (int j = 0; j < k; j++) {
        if (prob(i % np, j) < 0.0) {
          wrong_param = true;
          break;
        }
        p_tot += prob(i % np, j);
      }
      
      for (int j = 0; j < k; j++) {
        if (x(i % n, j) < 0.0 || !isInteger(x(i % n, j))) {
          wrong_x = true;
        } else {
          sum_x += x(i % n, j);
          prod_xfac += lfactorial(x(i % n, j));
          prod_pow_px += log(prob(i % np, j) / p_tot) * x(i % n, j);
        }
      }
      
    }
    
    if (wrong_param) {
      Rcpp::warning("NaNs produced");
      p[i] = NAN; 
    } else if (sum_x < 0.0 || sum_x != size[i % ns] || wrong_x) {
      p[i] = R_NegInf;
    } else {
      p[i] = n_fac - prod_xfac + prod_pow_px;
    }
  }
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
  return p;
}


// [[Rcpp::export]]
NumericMatrix cpp_rmnom(
    const int n,
    const NumericVector& size,
    const NumericMatrix& prob
  ) {
  
  int k = prob.ncol();
  int np = prob.nrow();
  int ns = size.length();
  bool wrong_param, missings;
  double p_tmp, size_left, sum_p, p_tot;
  
  NumericMatrix x(n, k);
  
  for (int i = 0; i < n; i++) {
    
    size_left = size[i % ns];
    sum_p = 1.0;
    p_tot = 0.0;
    wrong_param = false;
    missings = false;
    
    if (ISNAN(size[i % ns]))
      missings = true;
    
    for (int j = 0; j < k; j++) {
      if (ISNAN(prob(i % np, j)) || ISNAN(x(i % n, j))) {
        missings = true;
        break;
      }
    }
    
    if (missings) {
      for (int j = 0; j < k; j++)
        x(i, j) = NA_REAL;
      continue;
    } 
    
    // TODO:
    // sort prob(i,_) first?
    
    if (size[i % ns] < 0.0 || floor(size[i % ns]) != size[i % ns]) {
      wrong_param = true;
    } else {
      
      for (int j = 0; j < k; j++) {
        if (prob(i % np, j) < 0.0) {
          wrong_param = true;
          break;
        }
        p_tot += prob(i % np, j);
      }
      
      for (int j = 0; j < k-1; j++) {
        p_tmp = prob(i % np, j)/p_tot;
        x(i, j) = R::rbinom(size_left, p_tmp/sum_p);
        size_left -= x(i, j);
        sum_p -= p_tmp;
      }
      
      x(i, k-1) = size_left;
    }
    
    if (wrong_param) {
      Rcpp::warning("NaNs produced");
      for (int j = 0; j < k; j++)
        x(i, j) = NAN;
    }
  }
  
  return x;
}

