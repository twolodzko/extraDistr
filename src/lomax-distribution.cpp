#include <Rcpp.h>
using namespace Rcpp;


/*
*  Lomax distribution
*
*  Values:
*  x > 0
*
*  Parameters:
*  lambda > 0
*  kappa > 0
*
*  f(x)    = lambda*kappa / (1+lambda*x)^(kappa+1)
*  F(x)    = 1-(1+lambda*x)^-kappa
*  F^-1(p) = ((1-p)^(-1/kappa)-1) / lambda
*
*/

double pdf_lomax(double x, double lambda, double kappa) {
  if (lambda <= 0 || kappa <= 0)
    return NAN;
  if (x > 0)
    return lambda*kappa / std::pow(1+lambda*x, kappa+1);
  else
    return 0;
}

double cdf_lomax(double x, double lambda, double kappa) {
  if (lambda <= 0 || kappa <= 0)
    return NAN;
  if (x > 0)
    return 1-std::pow(1+lambda*x, -kappa);
  else
    return 0;
}

double invcdf_lomax(double p, double lambda, double kappa) {
  if (lambda <= 0 || kappa <= 0 || p < 0 || p > 1)
    return NAN;
  return (std::pow(1-p, -1/kappa)-1) / lambda;
}

double logpdf_lomax(double x, double lambda, double kappa) {
  if (lambda <= 0 || kappa <= 0)
    return NAN;
  if (x > 0)
    return std::log(lambda) + std::log(kappa) - std::log(1+lambda*x)*(kappa+1);
  else
    return -INFINITY;
}


// [[Rcpp::export]]
NumericVector cpp_dlomax(NumericVector x,
                         NumericVector lambda, NumericVector kappa,
                         bool log_prob = false) {

  int n = x.length();
  int nl = lambda.length();
  int nk = kappa.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nl, nk));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_lomax(x[i % n], lambda[i % nl], kappa[i % nk]);

  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::exp(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_plomax(NumericVector x,
                         NumericVector lambda, NumericVector kappa,
                         bool lower_tail = true, bool log_prob = false) {

  int n  = x.length();
  int nl = lambda.length();
  int nk = kappa.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nl, nk));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_lomax(x[i % n], lambda[i % nl], kappa[i % nk]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qlomax(NumericVector p,
                         NumericVector lambda, NumericVector kappa,
                         bool lower_tail = true, bool log_prob = false) {

  int n  = p.length();
  int nl = lambda.length();
  int nk = kappa.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nl, nk));
  NumericVector q(Nmax);

  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = std::exp(p[i]);

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      p[i] = 1-p[i];

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_lomax(p[i % n], lambda[i % nl], kappa[i % nk]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rlomax(int n,
                         NumericVector lambda, NumericVector kappa) {

  double u;
  int nl = lambda.length();
  int nk = kappa.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    u = R::runif(0, 1);
    x[i] = invcdf_lomax(u, lambda[i % nl], kappa[i % nk]);
  }

  return x;
}

