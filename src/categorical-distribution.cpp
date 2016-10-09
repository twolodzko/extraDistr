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
    bool log_prob = false
  ) {
  
  int n  = x.length();
  int np = prob.nrow();
  int Nmax = Rcpp::max(IntegerVector::create(n, np));
  int k = prob.ncol();
  NumericVector p(Nmax);
  bool missings;
  
  for (int i = 0; i < Nmax; i++) {
    
    missings = false;
    
    if (ISNAN(x[i]))
      missings = true;
    
    for (int j = 0; j < k; j++) {
      if (ISNAN(prob(i % np, j))) {
        missings = true;
        break;
      }
    }
    
    if (missings) {
      p[i] = NA_REAL;
      continue;
    }
    
    if (!isInteger(x[i]) || x[i] < 1.0 || x[i] > static_cast<double>(k)) {
      p[i] = 0.0;
    } else {
      double p_tot = 0.0;
      bool wrong_param = false;
      for (int j = 0; j < k; j++) {
        if (prob(i % np, j) < 0.0) {
          wrong_param = true;
          break;
        }
        p_tot += prob(i % np, j);
      }
      
      if (wrong_param) {
        Rcpp::warning("NaNs produced");
        p[i] = NAN;
      } else {
        p[i] = prob(i % np, static_cast<int>(x[i] - 1.0)) / p_tot;
      }
    }
  }

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
    
    return p;
}


// [[Rcpp::export]]
NumericVector cpp_pcat(
    const NumericVector& x,
    const NumericMatrix& prob,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = x.length();
  int np = prob.nrow();
  int Nmax = Rcpp::max(IntegerVector::create(n, np));
  int k = prob.ncol();
  NumericVector p(Nmax);
  bool missings;
  
  for (int i = 0; i < Nmax; i++) {
    
    missings = false;
    
    if (ISNAN(x[i]))
      missings = true;
    
    for (int j = 0; j < k; j++) {
      if (ISNAN(prob(i % np, j))) {
        missings = true;
        break;
      }
    }
    
    if (missings) {
      p[i] = NA_REAL;
      continue;
    }
    
    if (x[i] < 1.0) {
      p[i] = 0.0;
    } else if (x[i] > static_cast<double>(k)) {
      p[i] = 1.0;
    } else {
      bool wrong_param = false;
      p[i] = 0.0;
      int j = 0;
      while (j < static_cast<int>(x[i])) {
        if (prob(i % np, j) < 0.0) {
          wrong_param = true;
          break;
        }
        p[i] += prob(i % np, j);
        j++;
      }
      double p_tot = p[i];
      if (!wrong_param) {
        while (j < k) {
          if (prob(i % np, j) < 0.0) {
            wrong_param = true;
            break;
          }
          p_tot += prob(i % np, j);
          j++;
        }
      }

      if (wrong_param) {
        Rcpp::warning("NaNs produced");
        p[i] = NAN;
      } else {
        p[i] = p[i] / p_tot;
      }
      
    }
  }

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];
    
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
      
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qcat(
    const NumericVector& p,
    const NumericMatrix& prob,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = p.length();
  int np = prob.nrow();
  int Nmax = Rcpp::max(IntegerVector::create(n, np));
  int k = prob.ncol();
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);
    
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1.0 - pp[i];
  
  int jj;
  double pp_norm, p_tmp, p_tot;
  bool wrong_param, missings;
    
  for (int i = 0; i < Nmax; i++) {
    
    missings = false;
    
    if (ISNAN(pp[i]))
      missings = true;
    
    for (int j = 0; j < k; j++) {
      if (ISNAN(prob(i % np, j))) {
        missings = true;
        break;
      }
    }
    
    if (missings) {
      q[i] = NA_REAL;
      continue;
    }
    
    if (pp[i] < 0.0 || pp[i] > 1.0) {
      Rcpp::warning("NaNs produced");
      q[i] = NAN;
    } else if (pp[i] == 0.0) {
      q[i] = 1.0; 
    } else {
      pp_norm = pp[i];
      wrong_param = false;
      p_tmp = 1.0;
      p_tot = 0.0;
      jj = 0;
      
      for (int j = 0; j < k; j++) {
        if (prob(i % np, j) < 0.0) {
          wrong_param = true;
          break;
        }
        p_tot += prob(i % np, j);
      }
      
      for (int j = k-1; j >= 0; j--) {
        p_tmp -= prob(i % np, j) / p_tot;
        if (pp_norm > p_tmp) {
          jj = j;
          break;
        }
      }

      if (wrong_param) {
        Rcpp::warning("NaNs produced");
        q[i] = NAN;
      } else {
        q[i] = static_cast<double>(jj+1);
      }
    }
  }
      
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rcat(
    const int n,
    const NumericMatrix& prob
  ) {
  
  int np = prob.nrow();
  int k = prob.ncol();
  NumericVector x(n);
  
  int jj;
  double u, p_tmp, p_tot;
  bool wrong_param, missings;

  for (int i = 0; i < n; i++) {
    
    missings = false;

    for (int j = 0; j < k; j++) {
      if (ISNAN(prob(i % np, j))) {
        missings = true;
        break;
      }
    }
    
    if (missings) {
      x[i] = NA_REAL;
      continue;
    }
    
    u = rng_unif();
    wrong_param = false;
    p_tmp = 1.0;
    jj = 0;
    p_tot = 0.0;
    
    for (int j = 0; j < k; j++) {
      if (prob(i % np, j) < 0.0) {
        wrong_param = true;
        break;
      }
      p_tot += prob(i % np, j);
    }
    
    for (int j = k-1; j >= 0; j--) {
      p_tmp -= prob(i % np, j) / p_tot;
      if (u > p_tmp) {
        jj = j;
        break;
      }
    }

    if (wrong_param) {
      Rcpp::warning("NaNs produced");
      x[i] = NAN;
    } else {
      x[i] = static_cast<double>(jj+1);
    }
  }
  
  return x;
}

