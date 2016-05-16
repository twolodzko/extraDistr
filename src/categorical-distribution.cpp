#include <Rcpp.h>
#include "const.h"
#include "shared.h"
using namespace Rcpp;


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
NumericVector cpp_dcat(NumericVector x, NumericMatrix prob,
                       bool log_prob = false) {
  
  int n  = x.length();
  int np = prob.nrow();
  int Nmax = Rcpp::max(IntegerVector::create(n, np));
  int k = prob.ncol();
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++) {
    double p_tot = 0;
    bool wrong_p = false;
    for (int j = 0; j < k; j++) {
      if (prob(i % np, j) < 0 || prob(i % np, j) > 1) {
        wrong_p = true;
        break;
      }
      p_tot += prob(i % np, j)*P_NORM_CONST;
    }

    if (!tol_equal(p_tot/P_NORM_CONST, 1) || wrong_p) {
      Rcpp::warning("NaNs produced");
      p[i] = NAN;
    } else if (!isInteger(x[i]) || x[i] < 1 || x[i] > k) {
      p[i] = 0;
    } else {
      p[i] = prob(i % np, x[i]-1);
    }
  }

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
    
    return p;
}


// [[Rcpp::export]]
NumericVector cpp_pcat(NumericVector x, NumericMatrix prob,
                       bool lower_tail = true, bool log_prob = false) {
  
  int n  = x.length();
  int np = prob.nrow();
  int Nmax = Rcpp::max(IntegerVector::create(n, np));
  int k = prob.ncol();
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++) {
    if (x[i] < 1) {
      p[i] = 0;
    } else if (x[i] > k) {
      p[i] = 1;
    } else {
      bool wrong_p = false;
      p[i] = 0;
      int j = 0;
      while (j < std::min((int)x[i], k)) {
        if (prob(i % np, j) < 0 || prob(i % np, j) > 1) {
          wrong_p = true;
          break;
        }
        p[i] += prob(i % np, j)*P_NORM_CONST;
        j++;
      }
      double p_tot = p[i];
      while (j < k) {
        if (prob(i % np, j) < 0 || prob(i % np, j) > 1) {
          wrong_p = true;
          break;
        }
        p_tot += prob(i % np, j)*P_NORM_CONST;
        j++;
      }
      if (!tol_equal(p_tot/P_NORM_CONST, 1) || wrong_p) {
        Rcpp::warning("NaNs produced");
        p[i] = NAN;
      } else {
        p[i] = p[i]/P_NORM_CONST;
      }
      
    }
  }

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];
    
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
      
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qcat(NumericVector p, NumericMatrix prob,
                       bool lower_tail = true, bool log_prob = false) {
  
  int n  = p.length();
  int np = prob.nrow();
  int Nmax = Rcpp::max(IntegerVector::create(n, np));
  int k = prob.ncol();
  NumericVector q(Nmax);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = exp(p[i]);
    
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      p[i] = 1-p[i];
    
  for (int i = 0; i < Nmax; i++) {
    if (p[i] < 0 || p[i] > 1) {
      Rcpp::warning("NaNs produced");
      q[i] = NAN;
    } else {
      bool wrong_p = false;
      double cs_prob = 0;
      int j = 0;
      while (cs_prob < p[i]*P_NORM_CONST && j < k) {
        if (prob(i % np, j) < 0 || prob(i % np, j) > 1) {
          wrong_p = true;
          break;
        }
        cs_prob += prob(i % np, j)*P_NORM_CONST;
        j++;
      }
      if (p[i] == 0)
        q[i] = 1;
      else
        q[i] = j;
      
      while (j < k) {
        if (prob(i % np, j) < 0 || prob(i % np, j) > 1) {
          wrong_p = true;
          break;
        }
        cs_prob += prob(i % np, j)*P_NORM_CONST;
        j++;
      } 
      if (!tol_equal(cs_prob/P_NORM_CONST, 1) || wrong_p) {
        Rcpp::warning("NaNs produced");
        q[i] = NAN;
      }
    }
  }
      
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rcat(int n, NumericMatrix prob) {
  
  double u;
  int np = prob.nrow();
  int k = prob.ncol();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    bool wrong_p = false;
    u = R::runif(0, P_NORM_CONST);

    double cs_prob = 0;
    int j = 0;
    while (cs_prob < u && j < k) {
      if (prob(i % np, j) < 0 || prob(i % np, j) > 1) {
        wrong_p = true;
        break;
      }
      cs_prob += prob(i % np, j)*P_NORM_CONST;
      j++;
    }
    x[i] = j;
    
    while (j < k) {
      if (prob(i % np, j) < 0 || prob(i % np, j) > 1) {
        wrong_p = true;
        break;
      }
      cs_prob += prob(i % np, j)*P_NORM_CONST;
      j++;
    } 
    
    if (!tol_equal(cs_prob/P_NORM_CONST, 1) || wrong_p) {
      Rcpp::warning("NaNs produced");
      x[i] = NAN;
    }
  }
  
  return x;
}

