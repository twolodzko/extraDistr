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
  
  std::vector<int> dims;
  dims.push_back(x.nrow());
  dims.push_back(n.nrow());
  dims.push_back(k.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
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
      if (!isInteger(n(i % dims[1], j), false) || n(i % dims[1], j) < 0.0)
        wrong_n = true;
      sum_x += x(i % dims[0], j);
      n_tot += n(i % dims[1], j);
    }
    
    if (ISNAN(sum_x + n_tot + k[i % dims[2]])) {
      p[i] = sum_x + n_tot + k[i % dims[2]];
      continue;
    } 
    
    if (wrong_n || k[i % dims[2]] < 0.0 || k[i % dims[2]] > n_tot ||
        !isInteger(k[i % dims[2]], false)) {
      throw_warning = true;
      p[i] = NAN;
      continue;
    }
    
    for (int j = 0; j < m; j++) {
      if (x(i % dims[0], j) > n(i % dims[1], j) || x(i % dims[0], j) < 0.0 ||
          !isInteger(x(i % dims[0], j))) {
        wrong_x = true;
      } else {
        lncx_prod += R::lchoose(n(i % dims[1], j), x(i % dims[0], j));
      }
    }
    
    if (wrong_x || sum_x != k[i % dims[2]]) {
      p[i] = R_NegInf;
    } else {
      lNck = R::lchoose(n_tot, k[i % dims[2]]);
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
  
  std::vector<int> dims;
  dims.push_back(n.nrow());
  dims.push_back(k.length());
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
      if (!isInteger(n(i % dims[0], j), false) ||
          n(i % dims[0], j) < 0.0 || ISNAN(n(i % dims[0], j))) {
        wrong_values = true;
        break;
      }
      n_otr[0] += n(i % dims[0], j);
    }
    
    if (wrong_values || ISNAN(k[i % dims[1]]) || ISNAN(n(i % dims[0], 0)) ||
        !isInteger(n(i % dims[0], 0), false) || n(i % dims[0], 0) < 0 ||
        (n_otr[0] + n(i % dims[0], 0)) < k[i % dims[1]] ||
        !isInteger(k[i % dims[1]], false) || k[i % dims[1]] < 0.0) {
      throw_warning = true;
      for (int j = 0; j < m; j++)
        x(i, j) = NA_REAL;
      continue;
    }
    
    for (int j = 1; j < m; j++)
      n_otr[j] = n_otr[j-1] - n(i % dims[0], j);
    
    k_left = k[i % dims[1]];
    x(i, 0) = R::rhyper(n(i % dims[0], 0), n_otr[0], k_left);
    k_left -= x(i, 0);
    
    if (m > 2) {
      for (int j = 1; j < m-1; j++) {
        x(i, j) = R::rhyper(n(i % dims[0], j), n_otr[j], k_left);
        k_left -= x(i, j);
      }
    }
    
    x(i, m-1) = k_left;
    
  }
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

