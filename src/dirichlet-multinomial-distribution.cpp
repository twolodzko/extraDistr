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
*  Dirichlet-multinomial (multivariate Polya) distribution
*
*  Values:
*  x > 0
*
*  Parameters:
*  n > 0
*  alpha > 0    (R^k where k >= 2)
*  
*  where:
*  sum(x) == n
*
*/


// [[Rcpp::export]]
NumericVector cpp_ddirmnom(
    const NumericMatrix& x,
    const NumericVector& size,
    const NumericMatrix& alpha,
    bool log_prob = false
  ) {
  
  int n = x.nrow();
  int m = x.ncol();
  int k = alpha.ncol();
  int na = alpha.nrow();
  int ns = size.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, ns));
  k = std::min(m, k);
  NumericVector p(Nmax);
  
  if (k < 2)
    Rcpp::stop("Number of columns in 'alpha' should be >= 2.");
  if (m != k)
    Rcpp::stop("Number of columns in 'x' does not equal number of columns in 'alpha'.");
  
  double prod_tmp, sum_alpha, sum_x;
  bool wrong_x, wrong_param, missings;
  
  for (int i = 0; i < Nmax; i++) {
    prod_tmp = 0.0;
    sum_alpha = 0.0;
    sum_x = 0.0;
    wrong_x = false;
    wrong_param = false;
    missings = false;
    
    if (ISNAN(size[i % ns]))
      missings = true;
    
    for (int j = 0; j < k; j++) {
      if (ISNAN(alpha(i % na, j)) || ISNAN(x(i % n, j))) {
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
        if (alpha(i % na, j) <= 0.0) {
          wrong_param = true;
          break;
        }
        if (x(i % n, j) < 0.0 || !isInteger(x(i % n, j))) {
          wrong_x = true;
          break;
        }
        sum_x += x(i % n, j);
        prod_tmp += R::lgammafn(x(i % n, j) + alpha(i % na, j)) -
          (lfactorial(x(i % n, j)) + R::lgammafn(alpha(i % na, j)));
        sum_alpha += alpha(i % na, j);
      }
    }
    
    if (wrong_param) {
      Rcpp::warning("NaNs produced");
      p[i] = NAN;
    } else if (sum_x < 0.0 || sum_x != size[i % ns] || wrong_x) {
      p[i] = -INFINITY;
    } else {
      p[i] = (lfactorial(size[i % ns]) + R::lgammafn(sum_alpha)) -
        R::lgammafn(size[i % ns] + sum_alpha) + prod_tmp;
    }
  }
  
  if (!log_prob)
    for (int i = 0; i < n; i++)
      p[i] = exp(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericMatrix cpp_rdirmnom(
    const int n,
    const NumericVector& size,
    const NumericMatrix& alpha
  ) {
  
  int k = alpha.ncol();
  int na = alpha.nrow();
  int ns = size.length();
  NumericMatrix x(n, k);
  
  if (k < 2)
    Rcpp::stop("Number of columns in 'alpha' should be >= 2.");
  
  double size_left, row_sum;
  bool wrong_param, missings;
  
  for (int i = 0; i < n; i++) {
    size_left = size[i % ns];
    row_sum = 0.0;
    wrong_param = false;
    missings = false;
    NumericVector pi(k);
    
    if (ISNAN(size[i % ns]))
      missings = true;
    
    for (int j = 0; j < k; j++) {
      if (ISNAN(alpha(i % na, j))) {
        missings = true;
        break;
      }
    }
    
    if (missings) {
      for (int j = 0; j < k; j++)
        x(i, j) = NA_REAL;
      continue;
    } 
    
    if (size[i % ns] < 0.0 || floor(size[i % ns]) != size[i % ns]) {
      wrong_param = true;
    } else {
      for (int j = 0; j < k; j++) {
        if (ISNAN(alpha(i % na, j))) {
          missings = true;
          break;
        }
        if (alpha(i % na, j) <= 0.0) {
          wrong_param = true;
          break;
        }
        pi[j] = R::rgamma(alpha(i % na, j), 1.0);
        row_sum += pi[j];
      }
    }
    
    if (wrong_param) {
      Rcpp::warning("NaNs produced");
      for (int j = 0; j < k; j++)
        x(i, j) = NAN;
    } else if (size[i % ns] == 0.0) {
      for (int j = 0; j < k; j++)
        x(i, j) = 0.0;
    } else {
      double sum_p = 1.0;
      double p_tmp;
      for (int j = 0; j < k-1; j++) {
        p_tmp = pi[j] / row_sum;
        x(i, j) = R::rbinom(size_left, p_tmp/sum_p);
        size_left -= x(i, j);
        sum_p -= p_tmp;
      }
      x(i, k-1) = size_left;
    }
  }
  
  return x;
}

