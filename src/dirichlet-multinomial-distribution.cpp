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
    const bool& log_prob = false
  ) {
  
  if (std::min({static_cast<int>(x.nrow()),
                static_cast<int>(x.ncol()),
                static_cast<int>(size.length()),
                static_cast<int>(alpha.nrow()),
                static_cast<int>(alpha.ncol())}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    static_cast<int>(x.nrow()),
    static_cast<int>(size.length()),
    static_cast<int>(alpha.nrow())
  });

  int m = x.ncol();
  int k = alpha.ncol();
  k = std::min(m, k);
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  if (k < 2)
    Rcpp::stop("number of columns in alpha should be >= 2");
  if (m != k)
    Rcpp::stop("number of columns in x does not equal number of columns in alpha");
  
  double prod_tmp, sum_alpha, sum_x;
  bool wrong_x, wrong_param;
  
  for (int i = 0; i < Nmax; i++) {
    
    prod_tmp = 0.0;
    sum_alpha = 0.0;
    sum_x = 0.0;
    wrong_x = false;
    wrong_param = false;
    
    for (int j = 0; j < k; j++) {
      if (GETM(alpha, i, j) <= 0.0)
        wrong_param = true;
      if (GETM(x, i, j) < 0.0 || !isInteger(GETM(x, i, j)))
        wrong_x = true;
      
      sum_x += GETM(x, i, j);
      sum_alpha += GETM(alpha, i, j);
    }
    
#ifdef IEEE_754
    if (ISNAN(sum_x + sum_alpha + GETV(size, i))) {
      p[i] = sum_x + sum_alpha + GETV(size, i);
      continue;
    } 
#endif
    
    if (wrong_param || GETV(size, i) < 0.0 || !isInteger(GETV(size, i), false)) {
      throw_warning = true;
      p[i] = NAN;
      continue;
    }
    
    if (sum_x < 0.0 || sum_x != GETV(size, i) || wrong_x) {
      p[i] = R_NegInf;
    } else {
      
      for (int j = 0; j < k; j++) {
        prod_tmp += R::lgammafn(GETM(x, i, j) + GETM(alpha, i, j)) -
          (lfactorial(GETM(x, i, j)) + R::lgammafn(GETM(alpha, i, j)));
      }
      
      p[i] = (lfactorial(GETV(size, i)) + R::lgammafn(sum_alpha)) -
        R::lgammafn(GETV(size, i) + sum_alpha) + prod_tmp;
    }
  }
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericMatrix cpp_rdirmnom(
    const int& n,
    const NumericVector& size,
    const NumericMatrix& alpha
  ) {
  
  if (std::min({static_cast<int>(size.length()),
                static_cast<int>(alpha.nrow()),
                static_cast<int>(alpha.ncol())}) < 1) {
    Rcpp::warning("NAs produced");
    NumericMatrix out(n, alpha.ncol());
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }
  
  int k = alpha.ncol();
  NumericMatrix x(n, k);
  
  bool throw_warning = false;
  
  if (k < 2)
    Rcpp::stop("Number of columns in alpha should be >= 2");
  
  double size_left, row_sum, sum_p, p_tmp, sum_alpha;
  bool wrong_values;
  
  for (int i = 0; i < n; i++) {
    size_left = GETV(size, i);
    row_sum = 0.0;
    wrong_values = false;
    NumericVector pi(k);
    sum_alpha = 0.0;
    
    for (int j = 0; j < k; j++) {
      sum_alpha += GETM(alpha, i, j);
      
      if (GETM(alpha, i, j) <= 0.0) {
        wrong_values = true;
        break;
      }

      pi[j] = R::rgamma(GETM(alpha, i, j), 1.0);
      row_sum += pi[j];
    }
    
    if (wrong_values || ISNAN(sum_alpha + GETV(size, i)) ||
        GETV(size, i) < 0.0 || !isInteger(GETV(size, i), false)) {
      throw_warning = true;
      for (int j = 0; j < k; j++)
        x(i, j) = NA_REAL;
      continue;
    }
    
    if (GETV(size, i) == 0.0) {
      for (int j = 0; j < k; j++)
        x(i, j) = 0.0;
      continue;
    } 
    
    sum_p = 1.0;
    
    for (int j = 0; j < k-1; j++) {
      if ( size_left > 0.0 ) {
        p_tmp = pi[j] / row_sum;
        x(i, j) = R::rbinom(size_left, trunc_p(p_tmp/sum_p));
        size_left -= x(i, j);
        sum_p -= p_tmp;
      } else {
        break;
      }
    }
    
    x(i, k-1) = size_left;
    
  }
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

