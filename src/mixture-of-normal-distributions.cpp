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


// [[Rcpp::export]]
NumericVector cpp_dmixnorm(
    const NumericVector& x,
    const NumericMatrix& mu,
    const NumericMatrix& sigma,
    const NumericMatrix& alpha,
    const bool& log_prob = false
  ) {
  
  int Nmax = std::max({
    x.length(),
    static_cast<long int>(mu.nrow()),
    static_cast<long int>(sigma.nrow()),
    static_cast<long int>(alpha.nrow())
  });
  int k = alpha.ncol();
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  if (k != mu.ncol() || k != sigma.ncol())
    Rcpp::stop("sizes of mu, sigma, and alpha do not match");
  
  bool wrong_param;
  double alpha_tot, nans_sum;
  
  for (int i = 0; i < Nmax; i++) {
    wrong_param = false;
    alpha_tot = 0.0;
    nans_sum = 0.0;
    p[i] = 0.0;
    
    for (int j = 0; j < k; j++) {
      if (GETM(alpha, i, j) < 0.0 || GETM(sigma, i, j) <= 0.0)
        wrong_param = true;
      nans_sum += GETM(mu, i, j) + GETM(sigma, i, j);
      alpha_tot += GETM(alpha, i, j);
    }
    
    if (ISNAN(nans_sum + alpha_tot + GETV(x, i))) {
      p[i] = nans_sum + alpha_tot + GETV(x, i);
      continue;
    }
    
    if (wrong_param) {
      throw_warning = true;
      p[i] = NAN;
      continue;
    }
    
    for (int j = 0; j < k; j++) {
      p[i] += (GETM(alpha, i, j) / alpha_tot) *
        R::dnorm(GETV(x, i), GETM(mu, i, j), GETM(sigma, i, j), false);
    }
  }
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
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
  
  int Nmax = std::max({
    x.length(),
    static_cast<long int>(mu.nrow()),
    static_cast<long int>(sigma.nrow()),
    static_cast<long int>(alpha.nrow())
  });
  int k = alpha.ncol();
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  if (k != mu.ncol() || k != sigma.ncol())
    Rcpp::stop("sizes of mu, sigma, and alpha do not match");
  
  bool wrong_param;
  double alpha_tot, nans_sum;
  
  for (int i = 0; i < Nmax; i++) {
    wrong_param = false;
    alpha_tot = 0.0;
    nans_sum = 0.0;
    p[i] = 0.0;
    
    for (int j = 0; j < k; j++) {
      if (GETM(alpha, i, j) < 0.0 || GETM(sigma, i, j) < 0.0) {
        wrong_param = true;
        break;
      }
      nans_sum += GETM(mu, i, j) + GETM(sigma, i, j);
      alpha_tot += GETM(alpha, i, j);
    }
    
    if (ISNAN(nans_sum + alpha_tot + GETV(x, i))) {
      p[i] = nans_sum + alpha_tot + GETV(x, i);
      continue;
    }
    
    if (wrong_param) {
      throw_warning = true;
      p[i] = NAN;
      continue;
    }
    
    for (int j = 0; j < k; j++) {
      p[i] += (GETM(alpha, i, j) / alpha_tot) *
        R::pnorm(GETV(x, i), GETM(mu, i, j), GETM(sigma, i, j), true, false);
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
NumericVector cpp_rmixnorm(
    const int& n,
    const NumericMatrix& mu,
    const NumericMatrix& sigma,
    const NumericMatrix& alpha
  ) {
  
  int k = alpha.ncol();
  NumericVector x(n);
  
  bool throw_warning = false;
  
  if (k != mu.ncol() || k != sigma.ncol())
    Rcpp::stop("sizes of mu, sigma, and alpha do not match");
  
  int jj;
  bool wrong_param;
  double alpha_tot, nans_sum, u, p_tmp;
  NumericVector prob(k);
  
  for (int i = 0; i < n; i++) {
    jj = 0;
    wrong_param = false;
    u = rng_unif();
    p_tmp = 1.0;
    alpha_tot = 0.0;
    nans_sum = 0.0;
    
    for (int j = 0; j < k; j++) {
      if (GETM(alpha, i, j) < 0.0 || GETM(sigma, i, j) < 0.0) {
        wrong_param = true;
        break;
      }
      nans_sum += GETM(mu, i, j) + GETM(sigma, i, j);
      alpha_tot += GETM(alpha, i, j);
    }
    
    if (ISNAN(nans_sum + alpha_tot) || wrong_param) {
      throw_warning = true;
      x[i] = NA_REAL;
      continue;
    }
    
    for (int j = k-1; j >= 0; j--) {
      p_tmp -= GETM(alpha, i, j) / alpha_tot;
      if (u > p_tmp) {
        jj = j;
        break;
      }
    }

    x[i] = R::rnorm(GETM(mu, i, jj), GETM(sigma, i, jj)); 
  }
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

