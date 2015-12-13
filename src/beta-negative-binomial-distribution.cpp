#include <Rcpp.h>
using namespace Rcpp;


/*
*  Beta-negative binomial distribution
*
*  Values:
*  x
*
*  Parameters:
*  r > 0
*  alpha > 0
*  beta > 0
*
*  f(k) = gamma(r+k)/(k! gamma(r)) * beta(alpha+r, beta+k)/beta(alpha, beta)
*
*/

double pmf_bnbinom(double k, double r, double alpha, double beta) {
  if (std::floor(k) != k)
    return 0;
  if (alpha <= 0 || beta <= 0 || std::floor(r) != r)
    return NAN;
  return (R::gammafn(r+k) / (R::gammafn(k+1) * R::gammafn(r))) *
          R::beta(alpha+r, beta+k) / R::beta(alpha, beta);
}

double rng_bnbinom(double r, double alpha, double beta) {
  if (alpha <= 0 || beta <= 0 || std::floor(r) != r)
    return NAN;
  double prob = R::rbeta(alpha, beta);
  return R::rnbinom(r, prob);
}

double logpmf_bnbinom(double k, double r, double alpha, double beta) {
  if (std::floor(k) != k)
    return -INFINITY;
  if (alpha <= 0 || beta <= 0 || std::floor(r) != r)
    return NAN;
  return (R::lgammafn(r+k) - (R::lgammafn(k+1) + R::lgammafn(r))) +
         R::lbeta(alpha+r, beta+k) - R::lbeta(alpha, beta);
}


// [[Rcpp::export]]
NumericVector cpp_dbnbinom(NumericVector x,
                           NumericVector size, NumericVector alpha, NumericVector beta,
                           bool log_prob = false) {

  int n = x.length();
  int nn = size.length();
  int na = alpha.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = logpmf_bnbinom(x[i % n], size[i % nn], alpha[i % na], beta[i % nb]);

  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::exp(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pbnbinom(NumericVector x,
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
      p_tmp += std::exp(logpmf_bnbinom(j, size[i % nn], alpha[i % na], beta[i % nb]));
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
NumericVector cpp_rbnbinom(int n,
                           NumericVector size, NumericVector alpha, NumericVector beta) {

  int nn = size.length();
  int na = alpha.length();
  int nb = beta.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_bnbinom(size[i % nn], alpha[i % na], beta[i % nb]);

  return x;
}

