#include <Rcpp.h>
using namespace Rcpp;


/*
 *  Frechet distribution
 *
 *  Values:
 *  x > mu
 *
 *  Parameters:
 *  lambda > 0
 *  mu
 *  sigma > 0
 *
 *  z       = (x-mu)/sigma
 *  f(x)    = lambda/sigma * z^{-1-lambda} * exp(-z^-lambda)
 *  F(x)    = exp(-z^-lambda)
 *  F^-1(p) = mu + sigma * -log(p)^{-1/lambda}
 *
 */

double pdf_frechet(double x, double lambda, double mu, double sigma) {
  if (lambda <= 0 || sigma <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x > mu) {
    double z = (x-mu)/sigma;
    return lambda/sigma * pow(z, -1-lambda) * exp(-pow(z, -lambda));
  } else {
    return 0;
  }
}

double cdf_frechet(double x, double lambda, double mu, double sigma) {
  if (lambda <= 0 || sigma <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x > mu) {
    double z = (x-mu)/sigma;
    return exp(pow(-z, -lambda));
  } else {
    return 0;
  }
}

double invcdf_frechet(double p, double lambda, double mu, double sigma) {
  if (lambda <= 0 || sigma <= 0 || p < 0 || p > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return mu + sigma * pow(-log(p), -1/lambda);
}


// [[Rcpp::export]]
NumericVector cpp_dfrechet(NumericVector x,
                           NumericVector lambda, NumericVector mu, NumericVector sigma,
                           bool log_prob = false) {

  int n  = x.length();
  int nl = lambda.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nl, nm, ns));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_frechet(x[i % n], lambda[i % nl], mu[i % nm], sigma[i % ns]);

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pfrechet(NumericVector x,
                           NumericVector lambda, NumericVector mu, NumericVector sigma,
                           bool lower_tail = true, bool log_prob = false) {

  int n  = x.length();
  int nl = lambda.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nl, nm, ns));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_frechet(x[i % n], lambda[i % nl], mu[i % nm], sigma[i % ns]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qfrechet(NumericVector p,
                           NumericVector lambda, NumericVector mu, NumericVector sigma,
                           bool lower_tail = true, bool log_prob = false) {

  int n  = p.length();
  int nl = lambda.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nl, nm, ns));
  NumericVector q(Nmax);

  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = exp(p[i]);

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      p[i] = 1-p[i];

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_frechet(p[i % n], lambda[i % nl], mu[i % nm], sigma[i % ns]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rfrechet(int n,
                           NumericVector lambda, NumericVector mu, NumericVector sigma) {

  double u;
  int nl = lambda.length();
  int nm = mu.length();
  int ns = sigma.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    u = R::runif(0, 1);
    x[i] = invcdf_frechet(u, lambda[i % nl], mu[i % nm], sigma[i % ns]);
  }

  return x;
}

