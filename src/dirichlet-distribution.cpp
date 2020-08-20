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
NumericVector cpp_ddirichlet(
    const NumericMatrix& x,
    const NumericMatrix& alpha,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.nrow(), x.ncol(),
                alpha.nrow(), alpha.ncol()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.nrow(),
    alpha.nrow()
  });
  int m = x.ncol();
  int k = alpha.ncol();
  k = std::min(m, k);
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  if (k < 2)
    Rcpp::stop("number of columns in alpha should be >= 2");
  if (m != k)
    Rcpp::stop("number of columns in x does not equal number of columns in alpha");
  
  double prod_gamma, sum_alpha, p_tmp, beta_const, sum_x;
  bool wrong_alpha, wrong_x;

  for (int i = 0; i < Nmax; i++) {
    
    wrong_alpha = false;
    wrong_x = false;
    sum_alpha = 0.0;
    sum_x = 0.0;
    
    for (int j = 0; j < m; j++) {
      sum_alpha += GETM(alpha, i, j);
      sum_x += GETM(x, i, j);
      if (GETM(alpha, i, j) <= 0.0)
        wrong_alpha = true;
      if (GETM(x, i, j) < 0.0 || GETM(x, i, j) > 1.0)
        wrong_x = true;
    }
    
#ifdef IEEE_754
    if (ISNAN(sum_x + sum_alpha)) {
      p[i] = sum_x + sum_alpha;
      continue;
    }
#endif
    
    if (wrong_alpha) {
      throw_warning = true;
      p[i] = NAN;
    } else if (wrong_x) {
      p[i] = R_NegInf;
    } else {
      
      prod_gamma = 0.0;
      p_tmp = 0.0;
      
      for (int j = 0; j < m; j++) {
        prod_gamma += R::lgammafn(GETM(alpha, i, j));
        p_tmp += log(GETM(x, i, j)) * (GETM(alpha, i, j) - 1.0);
        
        if (GETM(alpha, i, j) == 1.0 && GETM(x, i, j) == 0.0)
          p_tmp = R_NegInf;
      }
      
      beta_const = prod_gamma - R::lgammafn(sum_alpha);
      p[i] = p_tmp - beta_const;
    }
  }

  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericMatrix cpp_rdirichlet(
    const int& n,
    const NumericMatrix& alpha
  ) {
  
  if (std::min({alpha.nrow(), alpha.ncol()}) < 1) {
    Rcpp::warning("NAs produced");
    NumericMatrix out(n, alpha.ncol());
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }

  int k = alpha.ncol();
  NumericMatrix x(n, k);
  
  bool throw_warning = false;
  
  if (k < 2)
    Rcpp::stop("number of columns in alpha should be >= 2");
  
  double row_sum, sum_alpha;
  bool wrong_values;

  for (int i = 0; i < n; i++) {
    sum_alpha = 0.0;
    row_sum = 0.0;
    wrong_values = false;

    for (int j = 0; j < k; j++) {
      sum_alpha += GETM(alpha, i, j);
      if (GETM(alpha, i, j) <= 0.0) {
        wrong_values = true;
        break;
      }
      
      x(i, j) = R::rgamma(GETM(alpha, i, j), 1.0);
      row_sum += x(i, j);
    }

    if (ISNAN(sum_alpha) || wrong_values) {
      throw_warning = true;
      for (int j = 0; j < k; j++)
        x(i, j) = NA_REAL;
    } else {
      for (int j = 0; j < k; j++)
        x(i, j) /= row_sum;
    }
  }
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

