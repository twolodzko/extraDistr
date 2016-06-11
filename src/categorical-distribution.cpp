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
  
  for (int i = 0; i < Nmax; i++) {
    double p_tot = 0.0;
    bool wrong_param = false;
    for (int j = 0; j < k; j++) {
      if (prob(i % np, j) < 0.0 || prob(i % np, j) > 1.0) {
        wrong_param = true;
        break;
      }
      p_tot += prob(i % np, j)*P_NORM_CONST;
    }

    if (!tol_equal(p_tot/P_NORM_CONST, 1.0) || wrong_param) {
      Rcpp::warning("NaNs produced");
      p[i] = NAN;
    } else if (!isInteger(x[i]) || x[i] < 1.0 || x[i] > k) {
      p[i] = 0.0;
    } else {
      p[i] = prob(i % np, x[i] - 1.0);
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
  
  for (int i = 0; i < Nmax; i++) {
    if (x[i] < 1.0) {
      p[i] = 0.0;
    } else if (static_cast<int>(x[i]) > k) {
      p[i] = 1.0;
    } else {
      bool wrong_param = false;
      p[i] = 0.0;
      int j = 0;
      while (j < static_cast<int>(x[i])) {
        if (prob(i % np, j) < 0.0 || prob(i % np, j) > 1.0) {
          wrong_param = true;
          break;
        }
        p[i] += prob(i % np, j)*P_NORM_CONST;
        j++;
      }
      double p_tot = p[i];
      if (!wrong_param) {
        while (j < k) {
          if (prob(i % np, j) < 0.0 || prob(i % np, j) > 1.0) {
            wrong_param = true;
            break;
          }
          p_tot += prob(i % np, j)*P_NORM_CONST;
          j++;
        }
      }

      if (wrong_param || !tol_equal(p_tot/P_NORM_CONST, 1.0)) {
        Rcpp::warning("NaNs produced");
        p[i] = NAN;
      } else {
        p[i] = p[i]/P_NORM_CONST;
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
  double pp_norm, p_tmp;
  bool wrong_param;
    
  for (int i = 0; i < Nmax; i++) {
    if (pp[i] < 0.0 || pp[i] > 1.0) {
      Rcpp::warning("NaNs produced");
      q[i] = NAN;
    } else if (pp[i] == 0.0) {
      q[i] = 1.0; 
    } else {
      pp_norm = pp[i]*P_NORM_CONST;
      wrong_param = false;
      p_tmp = P_NORM_CONST;
      jj = 0;
      
      for (int j = k-1; j >= 0; j--) {
        p_tmp -= prob(i % np, j)*P_NORM_CONST;
        if (prob(i % np, j) < 0.0 || prob(i % np, j) > 1.0) {
          wrong_param = true;
          break;
        }
        if (pp_norm > p_tmp) {
          jj = j;
          break;
        }
      }
      
      if (!wrong_param && jj > 0) {
        for (int j = jj-1; j >= 0; j--) {
          p_tmp -= prob(i % np, j)*P_NORM_CONST;
          if (prob(i % np, j) < 0.0 || prob(i % np, j) > 1.0) {
            wrong_param = true;
            break;
          }
        } 
      }
      
      if (wrong_param || !tol_equal(p_tmp/P_NORM_CONST, 0.0)) {
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
  double u, p_tmp;
  bool wrong_param;

  for (int i = 0; i < n; i++) {
    u = R::runif(0.0, P_NORM_CONST);
    wrong_param = false;
    p_tmp = P_NORM_CONST;
    jj = 0;
    
    for (int j = k-1; j >= 0; j--) {
      p_tmp -= prob(i % np, j)*P_NORM_CONST;
      if (prob(i % np, j) < 0.0 || prob(i % np, j) > 1.0 || p_tmp < 0) {
        wrong_param = true;
        break;
      }
      if (u > p_tmp) {
        jj = j;
        break;
      }
    }
    
    if (!wrong_param && jj > 0) {
      for (int j = jj-1; j >= 0; j--) {
        p_tmp -= prob(i % np, j)*P_NORM_CONST;
        if (prob(i % np, j) < 0.0 || prob(i % np, j) > 1.0) {
          wrong_param = true;
          break;
        }
      } 
    }
    
    if (wrong_param || !tol_equal(p_tmp/P_NORM_CONST, 0.0)) {
      Rcpp::warning("NaNs produced");
      x[i] = NAN;
    } else {
      x[i] = static_cast<double>(jj+1);
    }
  }
  
  return x;
}

