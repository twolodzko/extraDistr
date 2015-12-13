#include <Rcpp.h>
using namespace Rcpp;


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
NumericVector cpp_dmvhyper(NumericMatrix x,
                           NumericMatrix n, NumericVector k,
                           bool log_prob = false) {
  
  int nx = x.nrow();
  int nn = n.nrow();
  int m = n.ncol();
  int nk = k.length();
  int Nmax = Rcpp::max(IntegerVector::create(nx, nn, nk));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++) {
    
    int N = 0;
    for (int j = 0; j < m; j++)
      N += n(i % nn, j);
    
    if (k[i % nk] > N) {
      p[i] = NAN;
    } else {
      
      double lNck = R::lchoose(N, k[i % nk]);
      double row_sum = 0;
      double lncx_prod = 0;
      bool wrong_x = false;
      
      for (int j = 0; j < m; j++) {
        if (x(i % nx, j) > n(i % nn, j) || x(i % nx, j) < 0 ||
            std::floor(x(i % nx, j)) != x(i % nx, j)) {
          wrong_x = true;
          break;
        }
        lncx_prod += R::lchoose(n(i % nn, j), x(i % nx, j));
        row_sum += x(i % nx, j);
      }
      
      if (N < k[i % nk]) {
        p[i] = NAN;
      } else if (wrong_x || row_sum != k[i % nk]) {
        p[i] = -INFINITY;
      } else {
        p[i] = lncx_prod - lNck;
      }
    }
  }
  
  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::exp(p[i]);
  
  return p;
}


// [[Rcpp::export]]
IntegerMatrix cpp_rmvhyper(int nn, NumericMatrix n, NumericVector k) {
  
  int nr = n.nrow();
  int m = n.ncol();
  int nk = k.length();
  IntegerMatrix x(nn, m);
  IntegerVector n_otr(m);

  for (int i = 0; i < nn; i++) {
    
    n_otr[0] = 0;
    for (int j = 1; j < m; j++)
      n_otr[0] += n(i % nr, j);
    
    for (int j = 1; j < m; j++)
      n_otr[j] = n_otr[j-1] - n(i % nr, j);
    
    int k_left = k[i % nk];
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

