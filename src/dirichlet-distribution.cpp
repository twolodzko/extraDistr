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
using Rcpp::NumericMatrix;


/*
 *  Dirichlet distribution
 *
 *  Values:
 *  x
 *
 *  Parameters:
 *  alpha > 0    (R^k where k >= 2)
 *
 *  f(x) = Gamma(sum(alpha)) / prod(Gamma(alpha)) * prod_k x[k]^{k-1}
 *
 */


// [[Rcpp::export]]
NumericVector cpp_ddirichlet(
    const NumericMatrix& x,
    const NumericMatrix& alpha,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.nrow());
  dims.push_back(alpha.nrow());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  int m = x.ncol();
  int k = alpha.ncol();
  k = std::min(m, k);
  NumericVector p(Nmax);
  
  if (k < 2)
    Rcpp::stop("Number of columns in alpha should be >= 2");
  if (m != k)
    Rcpp::stop("Number of columns in x does not equal number of columns in alpha");
  
  double prod_gamma, sum_alpha, p_tmp, beta_const;
  bool wrong_alpha, missings, wrong_x;

  for (int i = 0; i < Nmax; i++) {
    
    wrong_alpha = false;
    missings = false;
    wrong_x = false;
    
    for (int j = 0; j < m; j++) {
      if (ISNAN(alpha(i % dims[1], j)) || ISNAN(x(i % dims[0], j))) {
        missings = true;
        break;
      }
      if (alpha(i % dims[1], j) <= 0.0)
        wrong_alpha = true;
      if (x(i % dims[0], j) < 0.0 || x(i % dims[0], j) > 1.0)
        wrong_x = true;
    }
    
    if (missings) {
      p[i] = NA_REAL;
    } else if (wrong_alpha) {
      Rcpp::warning("NaNs produced");
      p[i] = NAN;
    } else if (wrong_x) {
      p[i] = R_NegInf;
    } else {
      
      prod_gamma = 0.0;
      sum_alpha = 0.0;
      p_tmp = 0.0;
      
      for (int j = 0; j < m; j++) {
        prod_gamma += R::lgammafn(alpha(i % dims[1], j));
        sum_alpha += alpha(i % dims[1], j);
        p_tmp += log(x(i % dims[0], j)) * (alpha(i % dims[1], j) - 1.0);
        
        if (alpha(i % dims[1], j) == 1.0 && x(i % dims[0], j) == 0.0)
          p_tmp = R_NegInf;
      }
      
      beta_const = prod_gamma - R::lgammafn(sum_alpha);
      p[i] = p_tmp - beta_const;
    }
  }

  if (!log_prob)
    p = Rcpp::exp(p);

  return p;
}


// [[Rcpp::export]]
NumericMatrix cpp_rdirichlet(
    const int& n,
    const NumericMatrix& alpha
  ) {

  int k = alpha.ncol();
  int dims = alpha.nrow();
  NumericMatrix x(n, k);
  
  if (k < 2)
    Rcpp::stop("Number of columns in alpha should be >= 2");
  
  double row_sum;
  bool throw_warning;

  for (int i = 0; i < n; i++) {
    row_sum = 0.0;
    throw_warning = false;

    for (int j = 0; j < k; j++) {
      if (ISNAN(alpha(i % dims, j)) || alpha(i % dims, j) <= 0.0) {
        throw_warning = true;
        break;
      }
      
      x(i, j) = R::rgamma(alpha(i % dims, j), 1.0);
      row_sum += x(i, j);
    }

    if (throw_warning) {
      Rcpp::warning("NAs produced");
      for (int j = 0; j < k; j++)
        x(i, j) = NA_REAL;
    } else {
      for (int j = 0; j < k; j++)
        x(i, j) /= row_sum;
    }
  }

  return x;
}

