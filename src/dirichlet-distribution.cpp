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
                             NumericVector alpha,
                             bool log_prob = false) {

  if (is_true(any(alpha <= 0)))
    throw Rcpp::exception("Values of alpha be > 0.");
  if (alpha.length() < 2 && x.ncol() < 2)
    throw Rcpp::exception("Length of alpha and number of columns in x should be >= 2.");
  if (alpha.length() == x.ncol())
    throw Rcpp::exception("Length of alpha and number of columns in x should be the same.");

  int n = x.nrow();
  int k = alpha.length();
  NumericVector p(n);

  double prod_gamma = 0;
  double sum_alpha = 0;
  double p_tmp;

  for (int j = 0; j < k; j++) {
    prod_gamma += lgamma(alpha[j]);
    sum_alpha += alpha[j];
  }

  double beta_const = prod_gamma - lgamma(sum_alpha);

  for (int i = 0; i < n; i++) {
    if (x[i % n] >= 0 && x[i % n] <= 1) {
      p_tmp = 0;
      for (int j = 0; j < k; j++) {
        p_tmp += log(x(i, j)) * (alpha[j]-1);
        if (alpha[j] == 1 && x(i, j) == 0)
          p_tmp = -INFINITY;
      }
      p[i] = p_tmp - beta_const;
    } else {
      p[i] = 0;
    }
  }

  if (!log_prob)
    for (int i = 0; i < n; i++)
      p[i] = exp(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pdirichlet(NumericMatrix x,
                             NumericVector alpha,
                             bool lower_tail = true, bool log_prob = false,
                             int nsim = 10000) {

  if (is_true(any(alpha <= 0)))
    throw Rcpp::exception("Values of alpha be > 0.");
  if (alpha.length() < 2 && x.ncol() < 2)
    throw Rcpp::exception("Length of alpha and number of columns should be >= 2.");
  if (alpha.length() == x.ncol())
    throw Rcpp::exception("Length of alpha and number of columns should be the same.");

  double row_sum, p_tmp;
  int n = x.nrow();
  int k = alpha.length();
  NumericVector p(n, 0.0), y(k);

  for (int r = 0; r < nsim; r++) {
    for (int i = 0; i < n; i++) {
      row_sum = 0;
      p_tmp = 1.0;

      for (int j = 0; j < k; j++) {
        y[j] = R::rgamma(alpha[j], 1);
        row_sum += y[j];
      }

      for (int j = 0; j < k; j++) {
        y[j] = y[j] / row_sum;

        if (y[j] > x(i, j)) {
          p_tmp = 0;
          break;
        }
      }
      p[i] += p_tmp/nsim;
    }
  }

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericMatrix cpp_rdirichlet(int n,
                             NumericVector alpha) {

  if (is_true(any(alpha <= 0)))
    throw Rcpp::exception("Values of alpha be > 0.");
  if (alpha.length() < 2)
    throw Rcpp::exception("Length of alpha should be >= 2.");

  int k = alpha.length();
  NumericMatrix x(n, k);
  double row_sum;

  for (int i = 0; i < n; i++) {
    row_sum = 0;

    for (int j = 0; j < k; j++) {
      x(i, j) = R::rgamma(alpha[j], 1);
      row_sum += x(i, j);
    }

    for (int j = 0; j < k; j++)
      x(i, j) = x(i, j) / row_sum;
  }

  return x;
}

