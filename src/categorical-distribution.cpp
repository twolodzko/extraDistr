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
  bool missings, wrong_param;
  double p_tot;
  
  for (int i = 0; i < Nmax; i++) {
    
    p_tot = 0.0;
    wrong_param = false;
    missings = false;
    
    if (ISNAN(x[i % dims[0]]))
      missings = true;
    
    for (int j = 0; j < k; j++) {
      if (ISNAN(prob(i % dims[1], j)))
        missings = true;
      if (prob(i % dims[1], j) < 0.0)
        wrong_param = true;
      p_tot += prob(i % dims[1], j);
    }
    
    if (missings) {
      p[i] = NA_REAL;
      continue;
    }
    
    if (wrong_param) {
      Rcpp::warning("NaNs produced");
      p[i] = NAN;
      continue;
    }
    
    if (!isInteger(x[i % dims[0]]) || x[i % dims[0]] < 1.0 ||
        x[i % dims[0]] > static_cast<double>(k)) {
      p[i] = 0.0;
      continue;
    }
    
    p[i] = prob(i % dims[1], static_cast<int>(x[i % dims[0]] - 1.0)) / p_tot;
    
  }

  if (log_prob)
    p = Rcpp::log(p);
    
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
  bool missings, wrong_param;
  double p_tot;
  
  for (int i = 0; i < Nmax; i++) {
    
    p_tot = 0.0;
    wrong_param = false;
    missings = false;
    
    if (ISNAN(x[i % dims[0]]))
      missings = true;
    
    int j = 0;
    while (j < static_cast<int>(x[i % dims[0]])) {
      if (ISNAN(prob(i % dims[1], j)))
        missings = true;
      if (prob(i % dims[1], j) < 0.0)
        wrong_param = true;
      p_tot += prob(i % dims[1], j);
      j++;
    }
    
    p[i] = p_tot;
    
    if (!wrong_param) {
      while (j < k) {
        if (ISNAN(prob(i % dims[1], j)))
          missings = true;
        if (prob(i % dims[1], j) < 0.0)
          wrong_param = true;
        p_tot += prob(i % dims[1], j);
        j++;
      }
    }
    
    if (missings) {
      p[i] = NA_REAL;
      continue;
    }
    
    if (wrong_param) {
      Rcpp::warning("NaNs produced");
      p[i] = NAN;
      continue;
    }
    
    if (x[i] < 1.0) {
      p[i] = 0.0;
    } else if (x[i % dims[0]] >= static_cast<double>(k)) {
      p[i] = 1.0;
    } else {
      p[i] = p[i] / p_tot;
    }
    
  }

  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
      
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
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  int jj;
  double pp_norm, p_tmp, p_tot;
  bool wrong_param, missings;
    
  for (int i = 0; i < Nmax; i++) {
    
    p_tmp = 1.0;
    p_tot = 0.0;
    missings = false;
    wrong_param = false;
    
    if (ISNAN(pp[i % dims[0]]))
      missings = true;
    
    for (int j = 0; j < k; j++) {
      if (ISNAN(prob(i % dims[1], j)))
        missings = true;
      if (prob(i % dims[1], j) < 0.0)
        wrong_param = true;
      p_tot += prob(i % dims[1], j);
    }
    
    if (missings) {
      x[i] = NA_REAL;
      continue;
    }
    
    if (wrong_param || pp[i % dims[0]] < 0.0 || pp[i % dims[0]] > 1.0) {
      Rcpp::warning("NaNs produced");
      x[i] = NAN;
      continue;
    } 
    
    if (pp[i % dims[0]] == 0.0) {
      x[i] = 1.0; 
    } else {
      
      pp_norm = pp[i % dims[0]];
      jj = 0;
      
      for (int j = k-1; j >= 0; j--) {
        p_tmp -= prob(i % dims[1], j) / p_tot;
        if (pp_norm > p_tmp) {
          jj = j;
          break;
        }
      }

      x[i] = static_cast<double>(jj+1);
      
    }
  }
      
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
  double u, p_tmp, p_tot;
  bool wrong_param, missings;

  for (int i = 0; i < n; i++) {
    
    p_tot = 0.0;
    wrong_param = false;
    missings = false;

    for (int j = 0; j < k; j++) {
      if (ISNAN(prob(i % dims, j)))
        missings = true;
      if (prob(i % dims, j) < 0.0)
        wrong_param = true;
      p_tot += prob(i % dims, j);
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
    
    u = rng_unif();
    p_tmp = 1.0;
    jj = 0;
    
    for (int j = k-1; j >= 0; j--) {
      p_tmp -= prob(i % dims, j) / p_tot;
      if (u > p_tmp) {
        jj = j;
        break;
      }
    }
    
    x[i] = static_cast<double>(jj+1);
    
  }
  
  return x;
}

