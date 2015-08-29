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
  double z = std::abs(x-mu)/sigma;
  return 1/(2*sigma) * exp(-z);
}

double cdf_laplace(double x, double mu, double sigma) {
  double z = (x-mu)/sigma;
  if (x < mu)
    return exp(z)/2;
  else
    return 1 - exp(-z)/2;
}

double invcdf_laplace(double p, double mu, double sigma) {
  if (p < 0.5)
    return mu + sigma * log(2*p);
  else
    return mu - sigma * log(2*(1-p));
}

double rng_laplace(double mu, double sigma) {
  double u = R::runif(-0.5, 0.5);
  return mu + sigma * R::sign(u) * log(1 - 2*std::abs(u));
}


// [[Rcpp::export]]
NumericVector cpp_dlaplace(NumericVector x,
                           NumericVector mu, NumericVector sigma,
                           bool log_prob = false) {

  if (is_true(any(sigma <= 0)))
    throw Rcpp::exception("Values of sigma be > 0.");

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
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_plaplace(NumericVector x,
                           NumericVector mu, NumericVector sigma,
                           bool lower_tail = true, bool log_prob = false) {

  if (is_true(any(sigma <= 0)))
    throw Rcpp::exception("Values of sigma be > 0.");

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
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qlaplace(NumericVector p,
                           NumericVector mu, NumericVector sigma,
                           bool lower_tail = true, bool log_prob = false) {

  if (is_true(any(p < 0)) || is_true(any(p > 1)))
    throw Rcpp::exception("Probabilities should range from 0 to 1.");
  if (is_true(any(sigma <= 0)))
    throw Rcpp::exception("Values of sigma be > 0.");

  int n  = p.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns));
  NumericVector q(Nmax);

  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = exp(p[i]);

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

  if (is_true(any(sigma <= 0)))
    throw Rcpp::exception("Values of sigma be > 0.");

  int nm = mu.length();
  int ns = sigma.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_laplace(mu[i % nm], sigma[i % ns]);

  return x;
}

