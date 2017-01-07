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
NumericVector cpp_dmixnorm(
    const NumericVector& x,
    const NumericMatrix& mu,
    const NumericMatrix& sigma,
    const NumericMatrix& alpha,
    const bool& log_prob = false
) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mu.nrow());
  dims.push_back(sigma.nrow());
  dims.push_back(alpha.nrow());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  int k = alpha.ncol();
  NumericVector p(Nmax);
  
  if (k != mu.ncol() || k != sigma.ncol())
    Rcpp::stop("sizes of mu, sigma, and alpha do not match");
  
  bool wrong_param, missings;
  double alpha_tot;
  
  for (int i = 0; i < Nmax; i++) {
    wrong_param = false;
    alpha_tot = 0.0;
    p[i] = 0.0;
    missings = false;
    
    for (int j = 0; j < k; j++) {
      if (ISNAN(alpha(i % dims[3], j)) || ISNAN(mu(i % dims[1], j)) || ISNAN(sigma(i % dims[2], j))) {
        missings = true;
        break;
      }
      if (alpha(i % dims[3], j) < 0.0 || sigma(i % dims[2], j) <= 0.0) {
        wrong_param = true;
        break;
      }
      alpha_tot += alpha(i % dims[3], j);
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
      p[i] += (alpha(i % dims[3], j) / alpha_tot) * R::dnorm(x[i % dims[0]], mu(i % dims[1], j), sigma(i % dims[2], j), false);
  }
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pmixnorm(
    const NumericVector& x,
    const NumericMatrix& mu,
    const NumericMatrix& sigma,
    const NumericMatrix& alpha,
    const bool& lower_tail = true,
    const bool& log_prob = false
) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mu.nrow());
  dims.push_back(sigma.nrow());
  dims.push_back(alpha.nrow());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  int k = alpha.ncol();
  NumericVector p(Nmax);
  
  if (k != mu.ncol() || k != sigma.ncol())
    Rcpp::stop("sizes of mu, sigma, and alpha do not match");
  
  bool wrong_param, missings;
  double alpha_tot;
  
  for (int i = 0; i < Nmax; i++) {
    wrong_param = false;
    alpha_tot = 0.0;
    p[i] = 0.0;
    missings = false;
    
    for (int j = 0; j < k; j++) {
      if (ISNAN(alpha(i % dims[3], j)) || ISNAN(mu(i % dims[1], j)) || ISNAN(sigma(i % dims[2], j))) {
        missings = true;
        break;
      }
      if (alpha(i % dims[3], j) < 0.0 || sigma(i % dims[2], j) < 0.0) {
        wrong_param = true;
        break;
      }
      alpha_tot += alpha(i % dims[3], j);
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
      p[i] += (alpha(i % dims[3], j) / alpha_tot) * R::pnorm(x[i % dims[0]], mu(i % dims[1], j), sigma(i % dims[2], j), lower_tail, false);
  }
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rmixnorm(
    const int& n,
    const NumericMatrix& mu,
    const NumericMatrix& sigma,
    const NumericMatrix& alpha
) {
  
  std::vector<int> dims;
  dims.push_back(mu.nrow());
  dims.push_back(sigma.nrow());
  dims.push_back(alpha.nrow());
  int k = alpha.ncol();
  NumericVector x(n);
  
  if (k != mu.ncol() || k != sigma.ncol())
    Rcpp::stop("sizes of mu, sigma, and alpha do not match");
  
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
      if (ISNAN(alpha(i % dims[2], j)) || ISNAN(mu(i % dims[0], j)) || ISNAN(sigma(i % dims[1], j))) {
        missings = true;
        break;
      }
      if (alpha(i % dims[2], j) < 0.0 || sigma(i % dims[1], j) < 0.0) {
        wrong_param = true;
        break;
      }
      alpha_tot += alpha(i % dims[2], j);
    }
    
    if (missings || wrong_param) {
      Rcpp::warning("NAs produced");
      x[i] = NA_REAL;
      continue;
    }
    
    for (int j = k-1; j >= 0; j--) {
      p_tmp -= alpha(i % dims[2], j) / alpha_tot;
      if (u > p_tmp) {
        jj = j;
        break;
      }
    }

    x[i] = R::rnorm(mu(i % dims[0], jj), sigma(i % dims[1], jj)); 
  }
  
  return x;
}

