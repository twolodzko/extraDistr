#include <Rcpp.h>
using namespace Rcpp;


/*
*  Gamma-Poisson distribution
*
*  Values:
*  x >= 0
*
*  Parameters:
*  alpha > 0
*  beta > 0
*
*  f(x) = (gamma(x+beta)*alpha^x) / (gamma(beta)*(1+alpha)^(beta+x) * x!)
*
*/

double pmf_gpois(int x, double alpha, double beta) {
  if (x >= 0) {
    return (R::gammafn(x+beta)*pow(alpha, x)) /
           (R::gammafn(beta) * pow(1+alpha, beta+x) * R::gammafn(x+1));
  } else {
    return 0;
  }
}

double cdf_gpois(int x, double alpha, double beta) {
  double p_tmp = 0;
  for (int j = 0; j < x+1; j++)
    p_tmp += pmf_gpois(j, alpha, beta);
  return p_tmp;
}

double rng_gpois(double alpha, double beta) {
  double lambda = R::rgamma(alpha, beta);
  return R::rpois(lambda);
}

double logpmf_gpois(int x, double alpha, double beta) {
  if (x >= 0) {
    return (R::lgammafn(x+beta) + log(alpha)*x) -
           (R::lgammafn(beta) + log(1+alpha)*(beta+x) + R::lgammafn(x+1));
  } else {
    return -INFINITY;
  }
}

double cdf2_gpois(int x, double alpha, double beta) {
  double p_tmp = 0;
  for (int j = 0; j < x+1; j++)
    p_tmp += exp(logpmf_gpois(j, alpha, beta));
  return p_tmp;
}


// [[Rcpp::export]]
NumericVector cpp_dgpois(IntegerVector x,
                         NumericVector alpha, NumericVector beta,
                         bool log_prob = false) {

  if (is_true(any(alpha <= 0)) || is_true(any(beta <= 0)))
    throw Rcpp::exception("Values of alpha and beta should be > 0.");

  int n  = x.length();
  int na = alpha.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = logpmf_gpois(x[i % n], alpha[i % na], beta[i % nb]);

  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = exp(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pgpois(IntegerVector x,
                         NumericVector alpha, NumericVector beta,
                         bool lower_tail = true, bool log_prob = false) {

  if (is_true(any(alpha <= 0)) || is_true(any(beta <= 0)))
    throw Rcpp::exception("Values of alpha and beta should be > 0.");

  int n  = x.length();
  int na = alpha.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf2_gpois(x[i % n], alpha[i % na], beta[i % nb]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rgpois(int n,
                         NumericVector alpha, NumericVector beta) {

  if (is_true(any(alpha <= 0)) || is_true(any(beta <= 0)))
    throw Rcpp::exception("Values of alpha and beta should be > 0.");

  int na = alpha.length();
  int nb = beta.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_gpois(alpha[i % na], beta[i % nb]);

  return x;
}

