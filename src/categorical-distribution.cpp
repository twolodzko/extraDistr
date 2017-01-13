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
*  Categorical distribution
*
*  Values:
*  x
*
*  Parameters:
*  0 <= p <= 1
*  sum(p) = 1
*
*/


// [[Rcpp::export]]
NumericVector cpp_dcat(
    const NumericVector& x,
    const NumericMatrix& prob,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(prob.nrow());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  int k = prob.ncol();
  NumericVector p(Nmax);
  double p_tot;
  
  bool throw_warning = false;
  
  if (k < 2)
    Rcpp::stop("number of columns in prob is < 2");
  
  NumericMatrix prob_tab = Rcpp::clone(prob);
  
  for (int i = 0; i < dims[1]; i++) {
    p_tot = 0.0;
    for (int j = 0; j < k; j++) {
      p_tot += prob_tab(i, j);
      if (ISNAN(p_tot))
        break;
      if (prob_tab(i, j) < 0.0) {
        p_tot = NAN;
        throw_warning = true;
        break;
      }
    }
    for (int j = 0; j < k; j++)
      prob_tab(i, j) /= p_tot;
  }
  
  for (int i = 0; i < Nmax; i++) {
    if (ISNAN(GETV(x, i))) {
      p[i] = GETV(x, i);
      continue;
    }
    if (!isInteger(GETV(x, i)) || GETV(x, i) < 1.0 ||
        GETV(x, i) > TO_DBL(k)) {
      p[i] = 0.0;
      continue;
    }
    p[i] = prob_tab(i % dims[1], TO_INT(GETV(x, i)) - 1);
  }

  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
    
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pcat(
    const NumericVector& x,
    const NumericMatrix& prob,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(prob.nrow());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  int k = prob.ncol();
  NumericVector p(Nmax);
  double p_tot;
  
  bool throw_warning = false;
  
  if (k < 2)
    Rcpp::stop("number of columns in prob is < 2");
  
  NumericMatrix prob_tab = Rcpp::clone(prob);
  
  for (int i = 0; i < dims[1]; i++) {
    p_tot = 0.0;
    for (int j = 0; j < k; j++) {
      p_tot += prob_tab(i, j);
      if (ISNAN(p_tot))
        break;
      if (prob_tab(i, j) < 0.0) {
        p_tot = NAN;
        throw_warning = true;
        break;
      }
    }
    prob_tab(i, 0) /= p_tot;
    for (int j = 1; j < k; j++) {
      prob_tab(i, j) /= p_tot;
      prob_tab(i, j) += prob_tab(i, j-1);
    }
  }
  
  for (int i = 0; i < Nmax; i++) {
    if (ISNAN(GETV(x, i))) {
      p[i] = GETV(x, i);
      continue;
    }
    if (GETV(x, i) < 1.0) {
      p[i] = 0.0;
      continue;
    }
    if (GETV(x, i) >= TO_DBL(k)) {
      p[i] = 1.0;
      continue;
    }
    p[i] = prob_tab(i % dims[1], TO_INT(GETV(x, i)) - 1);
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
NumericVector cpp_qcat(
    const NumericVector& p,
    const NumericMatrix& prob,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(prob.nrow());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  int k = prob.ncol();
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  int jj;
  double p_tot;
  
  bool throw_warning = false;
  
  if (k < 2)
    Rcpp::stop("number of columns in prob is < 2");
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  NumericMatrix prob_tab = Rcpp::clone(prob);
  
  for (int i = 0; i < dims[1]; i++) {
    p_tot = 0.0;
    for (int j = 0; j < k; j++) {
      p_tot += prob_tab(i, j);
      if (ISNAN(p_tot))
        break;
      if (prob_tab(i, j) < 0.0) {
        p_tot = NAN;
        throw_warning = true;
        break;
      }
    }
    prob_tab(i, 0) /= p_tot;
    for (int j = 1; j < k; j++) {
      prob_tab(i, j) /= p_tot;
      prob_tab(i, j) += prob_tab(i, j-1);
    }
  }
  
  for (int i = 0; i < Nmax; i++) {
    if (ISNAN(GETV(p, i))) {
      x[i] = GETV(p, i);
      continue;
    }
    if (ISNAN(prob_tab(i % dims[1], 0))) {
      x[i] = prob_tab(i % dims[1], 0);
      continue;
    }
    if (GETV(p, i) < 0.0 || GETV(p, i) > 1.0) {
      x[i] = NAN;
      throw_warning = true;
      continue;
    }
    if (GETV(p, i) == 0.0) {
      x[i] = 1.0;
      continue;
    }
    if (GETV(p, i) == 1.0) {
      x[i] = TO_DBL(k);
      continue;
    }
    
    jj = 1;
    for (int j = 0; j < k; j++) {
      if (prob_tab(i % dims[1], j) >= GETV(p, i)) {
        jj = j+1;
        break;
      }
    }
    x[i] = TO_DBL(jj);
  }
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
      
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rcat(
    const int& n,
    const NumericMatrix& prob
  ) {
  
  int dims = prob.nrow();
  int k = prob.ncol();
  NumericVector x(n);
  int jj;
  double u, p_tot;
  
  bool throw_warning = false;
  
  if (k < 2)
    Rcpp::stop("number of columns in prob is < 2");

  NumericMatrix prob_tab = Rcpp::clone(prob);
  
  for (int i = 0; i < dims; i++) {
    p_tot = 0.0;
    for (int j = 0; j < k; j++) {
      p_tot += prob_tab(i, j);
      if (ISNAN(p_tot))
        break;
      if (prob_tab(i, j) < 0.0) {
        p_tot = NAN;
        throw_warning = true;
        break;
      }
    }
    prob_tab(i, 0) /= p_tot;
    for (int j = 1; j < k; j++) {
      prob_tab(i, j) /= p_tot;
      prob_tab(i, j) += prob_tab(i, j-1);
    }
  }
  
  for (int i = 0; i < n; i++) {
    if (ISNAN(prob_tab(i % dims, 0))) {
      x[i] = prob_tab(i % dims, 0);
      continue;
    }
    
    u = rng_unif();
    jj = 1;
    
    for (int j = 0; j < k; j++) {
      if (prob_tab(i % dims, j) >= u) {
        jj = j+1;
        break;
      }
    }
    x[i] = TO_DBL(jj);
  }
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

