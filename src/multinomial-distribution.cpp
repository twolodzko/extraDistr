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
* Multinomial distribution
* 
* x[i]       number of values of i-th category drawn
* n = sum(x) total number of draws
* p[i]       probability of drawing i-th category value
*  
* f(x) = n!/prod(x[i]!) * prod(p[i]^x[i])
*
*/


// [[Rcpp::export]]
NumericVector cpp_dmnom(
    const NumericMatrix& x,
    const NumericVector& size,
    const NumericMatrix& prob,
    const bool& log_prob = false
  ) {
  
  if (std::min({static_cast<int>(x.nrow()),
                static_cast<int>(x.ncol()),
                static_cast<int>(size.length()),
                static_cast<int>(prob.nrow()),
                static_cast<int>(prob.ncol())}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    static_cast<int>(x.nrow()),
    static_cast<int>(size.length()),
    static_cast<int>(prob.nrow())
  });
  int m = x.ncol();
  int k = prob.ncol();
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  if (m != k)
    Rcpp::stop("number of columns in x does not equal number of columns in prob");
  
  double n_fac, prod_xfac, prod_pow_px, sum_x, p_tot;
  bool wrong_param, wrong_x;
  
  for (int i = 0; i < Nmax; i++) {
    
    sum_x = 0.0;
    p_tot = 0.0;
    wrong_param = false;
    wrong_x = false;
    
    for (int j = 0; j < k; j++) {
      if (GETM(prob, i, j) < 0.0)
        wrong_param = true;
      p_tot += GETM(prob, i, j);
      sum_x += GETM(x, i, j);
    }
    
#ifdef IEEE_754
    if (ISNAN(p_tot + sum_x + GETV(size, i))) {
      p[i] = p_tot + sum_x + GETV(size, i);
      continue;
    } 
#endif
    
    if (wrong_param || GETV(size, i) < 0.0 ||
        !isInteger(GETV(size, i), false)) {
      throw_warning = true;
      p[i] = NAN; 
      continue;
    }
    
    prod_xfac = 0.0;
    prod_pow_px = 0.0;
    
    for (int j = 0; j < k; j++) {
      if (GETM(x, i, j) < 0.0 || !isInteger(GETM(x, i, j))) {
        wrong_x = true;
      } else {
        prod_xfac += lfactorial(GETM(x, i, j));
        prod_pow_px += log(GETM(prob, i, j) / p_tot) * GETM(x, i, j);
      }
    }
    
    if (wrong_x || sum_x < 0.0 || sum_x != GETV(size, i)) {
      p[i] = R_NegInf;
    } else {
      n_fac = lfactorial(GETV(size, i));
      p[i] = n_fac - prod_xfac + prod_pow_px;
    }
  }
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericMatrix cpp_rmnom(
    const int& n,
    const NumericVector& size,
    const NumericMatrix& prob
  ) {
  
  if (std::min({static_cast<int>(size.length()),
                static_cast<int>(prob.nrow()),
                static_cast<int>(prob.ncol())}) < 1) {
    Rcpp::warning("NAs produced");
    NumericMatrix out(n, prob.ncol());
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }
  
  int k = prob.ncol();
  bool wrong_values;
  double p_tmp, size_left, sum_p, p_tot;
  
  NumericMatrix x(n, k);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++) {
    
    size_left = GETV(size, i);
    sum_p = 1.0;
    p_tot = 0.0;
    wrong_values = false;
    
    // TODO:
    // sort prob(i,_) first?
    
    for (int j = 0; j < k; j++) {
      if (GETM(prob, i, j) < 0.0) {
        wrong_values = true;
        break;
      }
      p_tot += GETM(prob, i, j);
    }
    
    if (wrong_values || ISNAN(p_tot + GETV(size, i)) ||
        GETV(size, i) < 0.0 || !isInteger(GETV(size, i), false)) {
      throw_warning = true;
      for (int j = 0; j < k; j++)
        x(i, j) = NA_REAL;
      continue;
    }

    for (int j = 0; j < k-1; j++) {
      if ( size_left > 0.0 ) {
        p_tmp = GETM(prob, i, j)/p_tot;
        x(i, j) = R::rbinom(size_left, trunc_p(p_tmp/sum_p));
        size_left -= x(i, j);
        sum_p -= p_tmp;
      } else {
        break;
      }
    }
    
    x(i, k-1) = size_left;
    
  }
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

