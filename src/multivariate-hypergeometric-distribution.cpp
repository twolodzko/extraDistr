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
 * Multivariate hypergeometric distribution
 * 
 * n[i]       number of marbels from i-th category
 * N = sum(n) total number of marbels
 * x[i]       number of marbles of i-th category drawn
 * k = sum(x) number of marbles to be drawn
 * 
 * where k <= N and x[i] <= n[i]
 *  
 * f(x) = prod(choose(n, x)) / choose(N, k)
 *
 */


// [[Rcpp::export]]
NumericVector cpp_dmvhyper(
    const NumericMatrix& x,
    const NumericMatrix& n,
    const NumericVector& k,
    const bool& log_prob = false
  ) {
  
  if (std::min({static_cast<int>(x.nrow()),
                static_cast<int>(x.ncol()),
                static_cast<int>(n.nrow()),
                static_cast<int>(n.ncol()),
                static_cast<int>(k.length())}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    static_cast<int>(x.nrow()),
    static_cast<int>(n.nrow()),
    static_cast<int>(k.length())
  });
  int m = x.ncol();
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  if (m != n.ncol())
    Rcpp::stop("number of columns in x does not equal number of columns in n");

  bool wrong_n, wrong_x;
  double lNck, sum_x, lncx_prod, n_tot;
  
  for (int i = 0; i < Nmax; i++) {
    
    wrong_x = false;
    wrong_n = false;
    sum_x = 0.0;
    lncx_prod = 0.0;
    n_tot = 0.0;
    
    for (int j = 0; j < m; j++) {
      if (!isInteger(GETM(n, i, j), false) || GETM(n, i, j) < 0.0)
        wrong_n = true;
      sum_x += GETM(x, i, j);
      n_tot += GETM(n, i, j);
    }
    
#ifdef IEEE_754
    if (ISNAN(sum_x + n_tot + GETV(k, i))) {
      p[i] = sum_x + n_tot + GETV(k, i);
      continue;
    } 
#endif
    
    if (wrong_n || GETV(k, i) < 0.0 || GETV(k, i) > n_tot ||
        !isInteger(GETV(k, i), false)) {
      throw_warning = true;
      p[i] = NAN;
      continue;
    }
    
    for (int j = 0; j < m; j++) {
      if (GETM(x, i, j) > GETM(n, i, j) || GETM(x, i, j) < 0.0 ||
          !isInteger(GETM(x, i, j))) {
        wrong_x = true;
      } else {
        lncx_prod += R::lchoose(GETM(n, i, j), GETM(x, i, j));
      }
    }
    
    if (wrong_x || sum_x != GETV(k, i)) {
      p[i] = R_NegInf;
    } else {
      lNck = R::lchoose(n_tot, GETV(k, i));
      p[i] = lncx_prod - lNck;
    }
    
  }
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericMatrix cpp_rmvhyper(
    const int& nn,
    const NumericMatrix& n,
    const NumericVector& k
  ) {
  
  if (std::min({static_cast<int>(n.nrow()),
                static_cast<int>(n.ncol()),
                static_cast<int>(k.length())}) < 1) {
    Rcpp::warning("NAs produced");
    NumericMatrix out(nn, n.ncol());
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }
  
  int m = n.ncol();
  NumericMatrix x(nn, m);
  std::vector<double> n_otr(m);
  
  bool wrong_values;
  double k_left;
  
  bool throw_warning = false;

  for (int i = 0; i < nn; i++) {
    
    wrong_values = false;
    n_otr[0] = 0.0;
    
    for (int j = 1; j < m; j++) {
      if (!isInteger(GETM(n, i, j), false) ||
          GETM(n, i, j) < 0.0 || ISNAN(GETM(n, i, j))) {
        wrong_values = true;
        break;
      }
      n_otr[0] += GETM(n, i, j);
    }
    
    if (wrong_values || ISNAN(GETV(k, i)) || ISNAN(GETM(n, i, 0)) ||
        !isInteger(GETM(n, i, 0), false) || GETM(n, i, 0) < 0 ||
        (n_otr[0] + GETM(n, i, 0)) < GETV(k, i) ||
        !isInteger(GETV(k, i), false) || GETV(k, i) < 0.0) {
      throw_warning = true;
      for (int j = 0; j < m; j++)
        x(i, j) = NA_REAL;
      continue;
    }
    
    for (int j = 1; j < m; j++)
      n_otr[j] = n_otr[j-1] - GETM(n, i, j);
    
    k_left = GETV(k, i);
    x(i, 0) = R::rhyper(GETM(n, i, 0), n_otr[0], k_left);
    k_left -= x(i, 0);
    
    if (m > 2) {
      for (int j = 1; j < m-1; j++) {
        x(i, j) = R::rhyper(GETM(n, i, j), n_otr[j], k_left);
        k_left -= x(i, j);
      }
    }
    
    x(i, m-1) = k_left;
    
  }
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

