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
  
  if (std::min({x.length(), prob.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    static_cast<int>(x.length()),
    static_cast<int>(prob.nrow())
  });
  int k = prob.ncol();
  NumericVector p(Nmax);
  double p_tot;
  
  bool throw_warning = false;
  
  if (k < 2)
    Rcpp::stop("number of columns in prob is < 2");
  
  NumericMatrix prob_tab = Rcpp::clone(prob);
  
  for (int i = 0; i < prob.nrow(); i++) {
    p_tot = 0.0;
    for (int j = 0; j < k; j++) {
      p_tot += prob_tab(i, j);
#ifdef IEEE_754
      if (ISNAN(p_tot))
        break;
#endif
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
#ifdef IEEE_754
    if (ISNAN(GETV(x, i))) {
      p[i] = GETV(x, i);
      continue;
    }
#endif
    if (!isInteger(GETV(x, i)) || GETV(x, i) < 1.0 ||
        GETV(x, i) > to_dbl(k)) {
      p[i] = 0.0;
      continue;
    }
    if (is_large_int(GETV(x, i))) {
      Rcpp::warning("NAs introduced by coercion to integer range");
      p[i] = NA_REAL;
    }
    p[i] = GETM(prob_tab, i, to_pos_int(GETV(x, i)) - 1);
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
  
  if (std::min({x.length(), prob.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    static_cast<int>(x.length()),
    static_cast<int>(prob.nrow())
  });
  int k = prob.ncol();
  NumericVector p(Nmax);
  double p_tot;
  
  bool throw_warning = false;

  if (k < 2)
    Rcpp::stop("number of columns in prob is < 2");
  
  NumericMatrix prob_tab = Rcpp::clone(prob);
  
  for (int i = 0; i < prob.nrow(); i++) {
    p_tot = 0.0;
    for (int j = 0; j < k; j++) {
      p_tot += prob_tab(i, j);
#ifdef IEEE_754
      if (ISNAN(p_tot))
        break;
#endif
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
#ifdef IEEE_754
    if (ISNAN(GETV(x, i))) {
      p[i] = GETV(x, i);
      continue;
    }
#endif
    if (GETV(x, i) < 1.0) {
      p[i] = 0.0;
      continue;
    }
    if (GETV(x, i) >= to_dbl(k)) {
      p[i] = 1.0;
      continue;
    }
    if (is_large_int(GETV(x, i))) {
      Rcpp::warning("NAs introduced by coercion to integer range");
      p[i] = NA_REAL;
    }
    p[i] = GETM(prob_tab, i, to_pos_int(GETV(x, i)) - 1);
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
  
  if (std::min({p.length(), prob.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    static_cast<int>(p.length()),
    static_cast<int>(prob.nrow())
  });
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
  
  for (int i = 0; i < prob.nrow(); i++) {
    p_tot = 0.0;
    for (int j = 0; j < k; j++) {
      p_tot += prob_tab(i, j);
#ifdef IEEE_754
      if (ISNAN(p_tot))
        break;
#endif
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
#ifdef IEEE_754
    if (ISNAN(GETV(p, i))) {
      x[i] = GETV(p, i);
      continue;
    }
#endif
    if (ISNAN(GETM(prob_tab, i, 0))) {
      x[i] = GETM(prob_tab, i, 0);
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
      x[i] = to_dbl(k);
      continue;
    }
    
    jj = 1;
    for (int j = 0; j < k; j++) {
      if (GETM(prob_tab, i, j) >= GETV(p, i)) {
        jj = j+1;
        break;
      }
    }
    x[i] = to_dbl(jj);
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
  
  if (prob.length() < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  int k = prob.ncol();
  NumericVector x(n);
  int jj;
  double u, p_tot;
  
  bool throw_warning = false;
  
  if (k < 2)
    Rcpp::stop("number of columns in prob is < 2");

  NumericMatrix prob_tab = Rcpp::clone(prob);
  
  for (int i = 0; i < prob_tab.nrow(); i++) {
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
    if (ISNAN(GETM(prob_tab, i , 0))) {
      x[i] = GETM(prob_tab, i, 0);
      continue;
    }
    
    u = rng_unif();
    jj = 1;
    
    for (int j = 0; j < k; j++) {
      if (GETM(prob_tab, i, j) >= u) {
        jj = j+1;
        break;
      }
    }
    x[i] = to_dbl(jj);
  }
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

