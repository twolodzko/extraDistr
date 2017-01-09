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


// [[Rcpp::export]]
NumericVector cpp_dmixpois(
    const NumericVector& x,
    const NumericMatrix& lambda,
    const NumericMatrix& alpha,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(lambda.nrow());
  dims.push_back(alpha.nrow());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  int k = alpha.ncol();
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  if (k != lambda.ncol())
    Rcpp::stop("sizes of lambda and alpha do not match");
  
  bool wrong_param;
  double alpha_tot, nans_sum;
  
  for (int i = 0; i < Nmax; i++) {
    wrong_param = false;
    alpha_tot = 0.0;
    nans_sum = 0.0;
    p[i] = 0.0;
    
    for (int j = 0; j < k; j++) {
      if (alpha(i % dims[2], j) < 0.0 || lambda(i % dims[1], j) < 0.0) {
        wrong_param = true;
        break;
      }
      nans_sum += lambda(i % dims[1], j);
      alpha_tot += alpha(i % dims[2], j);
    }
    
    if (ISNAN(nans_sum + alpha_tot + x[i % dims[0]])) {
      p[i] = nans_sum + alpha_tot + x[i % dims[0]];
      continue;
    }
    
    if (wrong_param) {
      throw_warning = true;
      p[i] = NAN;
      continue;
    }
    
    if (x[i % dims[0]] < 0.0 || !isInteger(x[i % dims[0]])) {
      p[i] = 0.0;
      continue;
    }
    
    for (int j = 0; j < k; j++) {
      p[i] += (alpha(i % dims[2], j) / alpha_tot) *
        R::dpois(x[i % dims[0]], lambda(i % dims[1], j), false);
    }
  }
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pmixpois(
    const NumericVector& x,
    const NumericMatrix& lambda,
    const NumericMatrix& alpha,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(lambda.nrow());
  dims.push_back(alpha.nrow());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  int k = alpha.ncol();
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  if (k != lambda.ncol())
    Rcpp::stop("sizes of lambda and alpha do not match");
  
  bool wrong_param;
  double alpha_tot, nans_sum;
  
  for (int i = 0; i < Nmax; i++) {
    wrong_param = false;
    alpha_tot = 0.0;
    p[i] = 0.0;
    
    for (int j = 0; j < k; j++) {
      if (alpha(i % dims[2], j) < 0.0 || lambda(i % dims[1], j) < 0.0) {
        wrong_param = true;
        break;
      }
      nans_sum += lambda(i % dims[1], j);
      alpha_tot += alpha(i % dims[2], j);
    }
    
    if (ISNAN(nans_sum + alpha_tot + x[i % dims[0]])) {
      p[i] = nans_sum + alpha_tot + x[i % dims[0]];
      continue;
    }
    
    if (wrong_param) {
      throw_warning = true;
      p[i] = NAN;
      continue;
    }
    
    if (x[i % dims[0]] < 0.0) {
      p[i] = 0.0;
      continue;
    }
    
    for (int j = 0; j < k; j++) {
      p[i] += (alpha(i % dims[2], j) / alpha_tot) *
        R::ppois(x[i % dims[0]], lambda(i % dims[1], j), true, false);
    }
  }
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rmixpois(
    const int& n,
    const NumericMatrix& lambda,
    const NumericMatrix& alpha
  ) {
  
  std::vector<int> dims;
  dims.push_back(lambda.nrow());
  dims.push_back(alpha.nrow());
  int k = alpha.ncol();
  NumericVector x(n);
  
  bool throw_warning = false;
  
  if (k != lambda.ncol())
    Rcpp::stop("sizes of lambda and alpha do not match");
  
  int jj;
  bool wrong_param;
  double u, p_tmp, alpha_tot, nans_sum;
  NumericVector prob(k);
  
  for (int i = 0; i < n; i++) {
    jj = 0;
    wrong_param = false;
    u = rng_unif();
    p_tmp = 1.0;
    alpha_tot = 0.0;
    
    for (int j = 0; j < k; j++) {
      if (alpha(i % dims[1], j) < 0.0 || lambda(i % dims[0], j) < 0.0) {
        wrong_param = true;
        break;
      }
      nans_sum += lambda(i % dims[0], j);
      alpha_tot += alpha(i % dims[1], j);
    }
    
    if (ISNAN(nans_sum + alpha_tot) || wrong_param) {
      throw_warning = true;
      x[i] = NA_REAL;
      continue;
    }
    
    for (int j = k-1; j >= 0; j--) {
      p_tmp -= alpha(i % dims[1], j) / alpha_tot;
      if (u > p_tmp) {
        jj = j;
        break;
      }
    }
    
    x[i] = R::rpois(lambda(i % dims[0], jj)); 
  }
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

