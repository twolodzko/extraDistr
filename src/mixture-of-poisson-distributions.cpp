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


// [[Rcpp::export]]
NumericVector cpp_dmixpois(
    const NumericVector& x,
    const NumericMatrix& lambda,
    const NumericMatrix& alpha,
    const bool& log_prob = false
  ) {
  
  if (std::min({static_cast<int>(x.length()),
                static_cast<int>(lambda.nrow()),
                static_cast<int>(lambda.ncol()),
                static_cast<int>(alpha.nrow()),
                static_cast<int>(alpha.ncol())}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    static_cast<int>(x.length()),
    static_cast<int>(lambda.nrow()),
    static_cast<int>(alpha.nrow())
  });
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
      if (GETM(alpha, i, j) < 0.0 || GETM(lambda, i, j) < 0.0) {
        wrong_param = true;
        break;
      }
      nans_sum += GETM(lambda, i, j);
      alpha_tot += GETM(alpha, i, j);
    }
    
#ifdef IEEE_754
    if (ISNAN(nans_sum + alpha_tot + GETV(x, i))) {
      p[i] = nans_sum + alpha_tot + GETV(x, i);
      continue;
    }
#endif
    
    if (wrong_param) {
      throw_warning = true;
      p[i] = NAN;
      continue;
    }
    
    if (GETV(x, i) < 0.0 || !isInteger(GETV(x, i))) {
      p[i] = R_NegInf;
      continue;
    }
    
    if (!R_finite(GETV(x, i))) {
      p[i] = R_NegInf;
      continue;
    }
    
    double mx = R_NegInf;
    std::vector<double> tmp(k);
    
    for (int j = 0; j < k; j++) {
      // p[i] += (GETM(alpha, i, j) / alpha_tot) *
      //   R::dpois(GETV(x, i), GETM(lambda, i, j), false);
      tmp[j] = log(GETM(alpha, i, j)) - log(alpha_tot) +
        R::dpois(GETV(x, i), GETM(lambda, i, j), true);
      if (tmp[j] > mx)
        mx = tmp[j];
    }
    
    for (int j = 0; j < k; j++)
      p[i] += exp(tmp[j] - mx);  // log-sum-exp trick
    
    p[i] = log(p[i]) + mx;
    
  }
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
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
  
  if (std::min({static_cast<int>(x.length()),
                static_cast<int>(lambda.nrow()),
                static_cast<int>(lambda.ncol()),
                static_cast<int>(alpha.nrow()),
                static_cast<int>(alpha.ncol())}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    static_cast<int>(x.length()),
    static_cast<int>(lambda.nrow()),
    static_cast<int>(alpha.nrow())
  });
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
      if (GETM(alpha, i, j) < 0.0 || GETM(lambda, i, j) < 0.0) {
        wrong_param = true;
        break;
      }
      nans_sum += GETM(lambda, i, j);
      alpha_tot += GETM(alpha, i, j);
    }
    
#ifdef IEEE_754
    if (ISNAN(nans_sum + alpha_tot + GETV(x, i))) {
      p[i] = nans_sum + alpha_tot + GETV(x, i);
      continue;
    }
#endif
    
    if (wrong_param) {
      throw_warning = true;
      p[i] = NAN;
      continue;
    }
    
    if (GETV(x, i) < 0.0) {
      p[i] = lower_tail ? R_NegInf : 0.0; // here
      continue;
    }

    if (GETV(x, i) == R_PosInf) {
      p[i] = !lower_tail ? R_NegInf : 0.0;
      continue;
    }
    
    double mx = R_NegInf;
    std::vector<double> tmp(k);
    
    for (int j = 0; j < k; j++) {
      // p[i] += (GETM(alpha, i, j) / alpha_tot) *
      //   R::ppois(GETV(x, i), GETM(lambda, i, j), true, false);
      tmp[j] = log(GETM(alpha, i, j)) - log(alpha_tot) +
        R::ppois(GETV(x, i), GETM(lambda, i, j), lower_tail, true);
      if (tmp[j] > mx)
        mx = tmp[j];
    }
    
    for (int j = 0; j < k; j++)
      p[i] += exp(tmp[j] - mx);  // log-sum-exp trick
    
    p[i] = log(p[i]) + mx;
    
  }
  
  // if (!lower_tail)
  //   p = 1.0 - p;
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
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
  
  if (std::min({static_cast<int>(lambda.nrow()),
                static_cast<int>(lambda.ncol()),
                static_cast<int>(alpha.nrow()),
                static_cast<int>(alpha.ncol())}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
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
    nans_sum = 0.0;
    
    for (int j = 0; j < k; j++) {
      if (GETM(alpha, i, j) < 0.0 || GETM(lambda, i, j) < 0.0) {
        wrong_param = true;
        break;
      }
      nans_sum += GETM(lambda, i, j);
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
    
    x[i] = R::rpois(GETM(lambda, i, jj)); 
  }
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

