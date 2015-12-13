#include <Rcpp.h>
using namespace Rcpp;

/*
*  Beta-binomial distribution
*
*  Values:
*  x
*
*  Parameters:
*  k > 0
*  alpha > 0
*  beta > 0
*
*  f(k) = choose(n, k) * (beta(k+alpha, n-k+beta)) / (beta(alpha, beta))
*
*/

double pmf_bbinom(double k, double n, double alpha, double beta) {
  if (std::floor(k) != k)
    return 0;
  if (alpha < 0 || beta < 0 || std::floor(n) != n)
    return NAN;
  return R::choose(n, k) * R::beta(k+alpha, n-k+beta) / R::beta(alpha, beta);
}

double rng_bbinom(double n, double alpha, double beta) {
  if (alpha < 0 || beta < 0 || std::floor(n) != n)
    return NAN;
  double prob = R::rbeta(alpha, beta);
  return R::rbinom(n, prob);
}

double logpmf_bbinom(double k, double n, double alpha, double beta) {
  if (std::floor(k) != k)
    return -INFINITY;
  if (alpha < 0 || beta < 0 || std::floor(n) != n)
    return NAN;
  return R::lchoose(n, k) + R::lbeta(k+alpha, n-k+beta) - R::lbeta(alpha, beta);
}


// [[Rcpp::export]]
NumericVector cpp_dbbinom(NumericVector x,
                          NumericVector size, NumericVector alpha, NumericVector beta,
                          bool log_prob = false) {

  int n = x.length();
  int nn = size.length();
  int na = alpha.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = logpmf_bbinom(x[i % n], size[i % nn], alpha[i % na], beta[i % nb]);

  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::exp(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pbbinom(NumericVector x,
                          NumericVector size, NumericVector alpha, NumericVector beta,
                          bool lower_tail = true, bool log_prob = false) {

  int n = x.length();
  int nn = size.length();
  int na = alpha.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++) {
    double p_tmp = 0;
    for (int j = 0; j < x[i % n]+1; j++)
      p_tmp += std::exp(logpmf_bbinom(j, size[i % nn], alpha[i % na], beta[i % nb]));
    p[i] = p_tmp;
  }

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rbbinom(int n,
                          NumericVector size, NumericVector alpha, NumericVector beta) {

  int nn = size.length();
  int na = alpha.length();
  int nb = beta.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_bbinom(size[i % nn], alpha[i % na], beta[i % nb]);

  return x;
}

