#include <Rcpp.h>
using namespace Rcpp;


/*
 *  Dirichlet distribution
 *
 *  Values:
 *  x
 *
 *  Parameters:
 *  alpha > 0    (R^k where k >= 2)
 *
 *  f(x) = Gamma(sum(alpha)) / prod(Gamma(alpha)) * prod_k x[k]^{k-1}
 *
 */


// [[Rcpp::export]]
NumericVector cpp_ddirichlet(NumericMatrix x,
                             NumericMatrix alpha,
                             bool log_prob = false) {

  int n = x.nrow();
  int m = x.ncol();
  int k = alpha.ncol();
  int na = alpha.nrow();
  int Nmax = Rcpp::max(IntegerVector::create(n, na));
  k = std::min(m, k);
  NumericVector p(Nmax);
  
  if (k < 2)
    Rcpp::stop("Number of columns in 'alpha' should be >= 2.");
  if (m != k)
    Rcpp::stop("Number of columns in 'x' does not equal number of columns in 'alpha'.");

  for (int i = 0; i < Nmax; i++) {
    double prod_gamma = 0;
    double sum_alpha = 0;
    double p_tmp = 0;
    
    for (int j = 0; j < m; j++) {
      if (alpha(i % na, j) <= 0) {
        p[i] = NAN;
        break;
      }
      if (x(i % n, j) < 0 || x(i % n, j) > 1) {
        p[i] = -INFINITY;
        break;
      }
      
      prod_gamma += R::lgammafn(alpha(i % na, j));
      sum_alpha += alpha(i % na, j);
      p_tmp += std::log(x(i % n, j)) * (alpha(i % na, j)-1);
      
      if (alpha(i % na, j) == 1 && x(i % n, j) == 0)
        p_tmp = -INFINITY;
    }
    
    double beta_const = prod_gamma - R::lgammafn(sum_alpha);
    p[i] = p_tmp - beta_const;
  }

  if (!log_prob)
    for (int i = 0; i < n; i++)
      p[i] = std::exp(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericMatrix cpp_rdirichlet(int n,
                             NumericMatrix alpha) {

  int k = alpha.ncol();
  int na = alpha.nrow();
  NumericMatrix x(n, k);
  
  if (k < 2)
    Rcpp::stop("Number of columns in 'alpha' should be >= 2.");

  for (int i = 0; i < n; i++) {
    double row_sum = 0;
    bool wrong_alpha = false;

    for (int j = 0; j < k; j++) {
      if (alpha(i % na, j) <= 0) {
        wrong_alpha = true;
        break;
      }
      
      x(i, j) = R::rgamma(alpha(i % na, j), 1);
      row_sum += x(i, j);
    }

    if (wrong_alpha) {
      for (int j = 0; j < k; j++)
        x(i, j) = NAN;
    } else {
      for (int j = 0; j < k; j++)
        x(i, j) = x(i, j) / row_sum;
    }
  }

  return x;
}

