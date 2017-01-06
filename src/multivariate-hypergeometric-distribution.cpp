#include <Rcpp.h>
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
  
  if (m != n.ncol())
    Rcpp::stop("number of columns in x does not equal number of columns in n");

  bool missings, wrong_n, wrong_x;
  double lNck, row_sum, lncx_prod, n_tot;
  
  for (int i = 0; i < Nmax; i++) {
    
    missings = false;
    wrong_x = false;
    wrong_n = false;
    row_sum = 0.0;
    lncx_prod = 0.0;
    n_tot = 0.0;
    
    for (int j = 0; j < m; j++) {
      if (ISNAN(x(i % dims[0], j)) || ISNAN(n(i % dims[1], j)))
        missings = true;
      if (!isInteger(n(i % dims[1], j), false) || n(i % dims[1], j) < 0.0)
        wrong_n = true;
      n_tot += n(i % dims[1], j);
    }
    
    if (missings || ISNAN(k[i % dims[2]])) {
      p[i] = NA_REAL;
      continue;
    } 
    
    if (wrong_n || k[i % dims[2]] < 0.0 || k[i % dims[2]] > n_tot ||
        !isInteger(k[i % dims[2]], false)) {
      Rcpp::warning("NaNs produced");
      p[i] = NAN;
      continue;
    }
    
    lNck = R::lchoose(n_tot, k[i % dims[2]]);
    
    for (int j = 0; j < m; j++) {
      if (x(i % dims[0], j) > n(i % dims[1], j) || x(i % dims[0], j) < 0.0 ||
          !isInteger(x(i % dims[0], j))) {
        wrong_x = true;
      } else {
        lncx_prod += R::lchoose(n(i % dims[1], j), x(i % dims[0], j));
        row_sum += x(i % dims[0], j);
      }
    }
    
    if (wrong_x || row_sum != k[i % dims[2]]) {
      p[i] = R_NegInf;
    } else {
      p[i] = lncx_prod - lNck;
    }
    
  }
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
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
  
  bool throw_warning;
  double k_left;

  for (int i = 0; i < nn; i++) {
    
    throw_warning = false;
    n_otr[0] = 0.0;
    
    for (int j = 1; j < m; j++) {
      if (!isInteger(n(i % dims[0], j), false) ||
          n(i % dims[0], j) < 0.0 || ISNAN(n(i % dims[0], j))) {
        throw_warning = true;
        break;
      }
      n_otr[0] += n(i % dims[0], j);
    }
    
    if (throw_warning || ISNAN(k[i % dims[1]]) || ISNAN(n(i % dims[0], 0)) ||
        !isInteger(n(i % dims[0], 0), false) || n(i % dims[0], 0) < 0 ||
        (n_otr[0] + n(i % dims[0], 0)) < k[i % dims[1]] ||
        !isInteger(k[i % dims[1]], false) || k[i % dims[1]] < 0.0) {
      Rcpp::warning("NAs produced");
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

  return x;
}

