#include <Rcpp.h>
using namespace Rcpp;


/*
 *  Laplace distribution
 *
 *  Values:
 *  x
 *
 *  Parameters:
 *  mu
 *  sigma > 0
 *
 *  z = (x-mu)/sigma
 *  f(x)    = 1/(2*sigma) * exp(-|z|)
 *  F(x)    = { 1/2 * exp(z)                 if   x < mu
 *            { 1 - 1/2 * exp(z)             otherwise
 *  F^-1(p) = { mu + sigma * log(2*p)        if p <= 0.5
 *            { mu + sigma * log(2*(1-p))    otherwise
 *
 */

double pdf_laplace(double x, double mu, double sigma) {
  if (sigma <= 0)
    return NAN;
  double z = std::abs(x-mu)/sigma;
  return 1/(2*sigma) * std::exp(-z);
}

double cdf_laplace(double x, double mu, double sigma) {
  if (sigma <= 0)
    return NAN;
  double z = (x-mu)/sigma;
  if (x < mu)
    return std::exp(z)/2;
  else
    return 1 - std::exp(-z)/2;
}

double invcdf_laplace(double p, double mu, double sigma) {
  if (sigma <= 0 || p < 0 || p > 1)
    return NAN;
  if (p < 0.5)
    return mu + sigma * std::log(2*p);
  else
    return mu - sigma * std::log(2*(1-p));
}

double rng_laplace(double mu, double sigma) {
  if (sigma <= 0)
    return NAN;
  double u = R::runif(-0.5, 0.5);
  return mu + sigma * R::sign(u) * std::log(1 - 2*std::abs(u));
}


// [[Rcpp::export]]
NumericVector cpp_dlaplace(NumericVector x,
                           NumericVector mu, NumericVector sigma,
                           bool log_prob = false) {

  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns));
  NumericVector p(Nmax);
  double z;

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_laplace(x[i % n], mu[i % nm], sigma[i % ns]);

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_plaplace(NumericVector x,
                           NumericVector mu, NumericVector sigma,
                           bool lower_tail = true, bool log_prob = false) {

  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns));
  NumericVector p(Nmax);
  double z;

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_laplace(x[i % n], mu[i % nm], sigma[i % ns]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qlaplace(NumericVector p,
                           NumericVector mu, NumericVector sigma,
                           bool lower_tail = true, bool log_prob = false) {

  int n  = p.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns));
  NumericVector q(Nmax);

  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = std::exp(p[i]);

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      p[i] = 1-p[i];

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_laplace(p[i % n], mu[i % nm], sigma[i % ns]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rlaplace(int n,
                           NumericVector mu, NumericVector sigma) {

  int nm = mu.length();
  int ns = sigma.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_laplace(mu[i % nm], sigma[i % ns]);

  return x;
}

