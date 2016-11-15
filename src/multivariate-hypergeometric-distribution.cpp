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
  
  int nx = x.nrow();
  int nr = n.nrow();
  int m = x.ncol();
  int nk = k.length();
  int Nmax = Rcpp::max(IntegerVector::create(nx, nr, nk));
  NumericVector p(Nmax);
  
  if (m != n.ncol())
    Rcpp::stop("Number of columns in 'x' does not equal number of columns in 'n'.");

  bool missings, wrong_n, wrong_x;
  double lNck, row_sum, lncx_prod, N;
  
  for (int i = 0; i < Nmax; i++) {
    
    missings = false;
    wrong_x = false;
    wrong_n = false;
    row_sum = 0.0;
    lncx_prod = 0.0;
    N = 0.0;
    
    for (int j = 0; j < m; j++) {
      if (ISNAN(x(i % nx, j)) || ISNAN(n(i % nr, j)))
        missings = true;
      if (floor(n(i % nr, j)) != n(i % nr, j) || n(i % nr, j) < 0.0)
        wrong_n = true;
      N += n(i % nr, j);
    }
    
    if (missings || ISNAN(k[i % nk])) {
      p[i] = NA_REAL;
      continue;
    } 
    
    if (wrong_n || k[i % nk] < 0.0 || k[i % nk] > N ||
        floor(k[i % nk]) != k[i % nk]) {
      Rcpp::warning("NaNs produced");
      p[i] = NAN;
      continue;
    }
    
    lNck = R::lchoose(N, k[i % nk]);
    
    for (int j = 0; j < m; j++) {
      if (x(i % nx, j) > n(i % nr, j) || x(i % nx, j) < 0.0 ||
          !isInteger(x(i % nx, j))) {
          wrong_x = true;
      } else {
        lncx_prod += R::lchoose(n(i % nr, j), x(i % nx, j));
        row_sum += x(i % nx, j);
      }
    }
    
    if (wrong_x || row_sum != k[i % nk]) {
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
  
  int nr = n.nrow();
  int m = n.ncol();
  int nk = k.length();
  NumericMatrix x(nn, m);
  std::vector<double> n_otr(m);
  
  bool missings, wrong_n;
  double k_left;

  for (int i = 0; i < nn; i++) {
    
    missings = false;
    wrong_n = false;
    n_otr[0] = 0.0;
    
    for (int j = 1; j < m; j++) {
      if (floor(n(i % nr, j)) != n(i % nr, j) || n(i % nr, j) < 0.0)
        wrong_n = true;
      if (ISNAN(n(i % nr, j)))
        missings = true;
      n_otr[0] += n(i % nr, j);
    }
    
    if (missings || ISNAN(k[i % nk]) || ISNAN(n(i % nr, 0))) {
      for (int j = 0; j < m; j++)
        x(i, j) = NA_REAL;
      continue;
    }
    
    if (wrong_n || floor(n(i % nr, 0)) != n(i % nr, 0) ||
        n(i % nr, 0) < 0 || (n_otr[0] + n(i % nr, 0)) < k[i % nk] ||
        floor(k[i % nk]) != k[i % nk] || k[i % nk] < 0.0) {
      Rcpp::warning("NaNs produced");
      for (int j = 0; j < m; j++)
        x(i, j) = NAN;
      continue;
    }
    
    for (int j = 1; j < m; j++)
      n_otr[j] = n_otr[j-1] - n(i % nr, j);
    
    k_left = k[i % nk];
    x(i, 0) = R::rhyper(n(i % nr, 0), n_otr[0], k_left);
    k_left -= x(i, 0);
    
    if (m > 2) {
      for (int j = 1; j < m-1; j++) {
        x(i, j) = R::rhyper(n(i % nr, j), n_otr[j], k_left);
        k_left -= x(i, j);
      }
    }
    
    x(i, m-1) = k_left;
    
  }

  return x;
}

