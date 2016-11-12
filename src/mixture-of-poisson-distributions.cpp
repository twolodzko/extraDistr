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


// [[Rcpp::export]]
NumericVector cpp_dmixpois(
    const NumericVector& x,
    const NumericMatrix& lambda,
    const NumericMatrix& alpha,
    bool log_prob = false
) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(lambda.nrow());
  dims.push_back(alpha.nrow());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  int k = alpha.ncol();
  NumericVector p(Nmax);
  
  if (k != lambda.ncol())
    Rcpp::stop("sizes of lambda and alpha do not match");
  
  bool wrong_param, missings;
  double alpha_tot;
  
  for (int i = 0; i < Nmax; i++) {
    wrong_param = false;
    alpha_tot = 0.0;
    p[i] = 0.0;
    missings = false;
    
    for (int j = 0; j < k; j++) {
      if (ISNAN(alpha(i % dims[2], j)) || ISNAN(lambda(i % dims[1], j))) {
        missings = true;
        break;
      }
      if (alpha(i % dims[2], j) < 0.0 || lambda(i % dims[1], j) < 0.0) {
        wrong_param = true;
        break;
      }
      alpha_tot += alpha(i % dims[2], j);
    }
    
    if (missings || ISNAN(x[i % dims[0]])) {
      p[i] = NA_REAL;
      continue;
    }
    
    if (wrong_param) {
      Rcpp::warning("NaNs produced");
      p[i] = NAN;
      continue;
    }
    
    for (int j = 0; j < k; j++)
      p[i] += (alpha(i % dims[2], j) / alpha_tot) * R::dpois(x[i % dims[0]], lambda(i % dims[1], j), false);
  }
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pmixpois(
    const NumericVector& x,
    const NumericMatrix& lambda,
    const NumericMatrix& alpha,
    bool lower_tail = true, bool log_prob = false
) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(lambda.nrow());
  dims.push_back(alpha.nrow());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  int k = alpha.ncol();
  NumericVector p(Nmax);
  
  if (k != lambda.ncol())
    Rcpp::stop("sizes of 'lambda' and 'alpha' do not match");
  
  bool wrong_param, missings;
  double alpha_tot;
  
  for (int i = 0; i < Nmax; i++) {
    wrong_param = false;
    alpha_tot = 0.0;
    p[i] = 0.0;
    missings = false;
    
    for (int j = 0; j < k; j++) {
      if (ISNAN(alpha(i % dims[2], j)) || ISNAN(lambda(i % dims[1], j))) {
        missings = true;
        break;
      }
      if (alpha(i % dims[2], j) < 0.0 || lambda(i % dims[1], j) < 0.0) {
        wrong_param = true;
        break;
      }
      alpha_tot += alpha(i % dims[2], j);
    }
    
    if (missings || ISNAN(x[i % dims[0]])) {
      p[i] = NA_REAL;
      continue;
    }
    
    if (wrong_param) {
      Rcpp::warning("NaNs produced");
      p[i] = NAN;
      continue;
    }
    
    for (int j = 0; j < k; j++)
      p[i] += (alpha(i % dims[2], j) / alpha_tot) * R::ppois(x[i % dims[0]], lambda(i % dims[1], j), lower_tail, false);
  }
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rmixpois(
    const int n,
    const NumericMatrix& lambda,
    const NumericMatrix& alpha
) {
  
  std::vector<int> dims;
  dims.push_back(lambda.nrow());
  dims.push_back(alpha.nrow());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  int k = alpha.ncol();
  NumericVector x(n);
  
  if (k != lambda.ncol())
    Rcpp::stop("sizes of lambda and alpha do not match");
  
  int jj;
  bool wrong_param, missings;
  double u, p_tmp, alpha_tot;
  NumericVector prob(k);
  
  for (int i = 0; i < n; i++) {
    jj = 0;
    wrong_param = false;
    u = rng_unif();
    p_tmp = 1.0;
    alpha_tot = 0.0;
    missings = false;
    
    for (int j = 0; j < k; j++) {
      if (ISNAN(alpha(i % dims[1], j)) || ISNAN(lambda(i % dims[0], j))) {
        missings = true;
        break;
      }
      if (alpha(i % dims[1], j) < 0.0 || lambda(i % dims[0], j) < 0.0) {
        wrong_param = true;
        break;
      }
      alpha_tot += alpha(i % dims[1], j);
    }
    
    if (missings) {
      x[i] = NA_REAL;
      continue;
    }
    
    if (wrong_param) {
      Rcpp::warning("NaNs produced");
      x[i] = NAN;
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
  
  return x;
}

