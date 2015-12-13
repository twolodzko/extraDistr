#include <Rcpp.h>
using namespace Rcpp;


/*
 *  Gumbel distribution
 *
 *  Values:
 *  x
 *
 *  Parameters:
 *  mu
 *  sigma > 0
 *
 *  z       = (x-mu)/sigma
 *  f(x)    = 1/sigma * exp(-(z+exp(-z)))
 *  F(x)    = exp(-exp(-z))
 *  F^-1(p) = mu - sigma * log(-log(p))
 *
 */

double pdf_gumbel(double x, double mu, double sigma) {
  if (sigma <= 0)
    return NAN;
  double z = (x-mu)/sigma;
  return 1/sigma * std::exp(-(z+std::exp(-z)));
}

double cdf_gumbel(double x, double mu, double sigma) {
  if (sigma <= 0)
    return NAN;
  double z = (x-mu)/sigma;
  return std::exp(-std::exp(-z));
}

double invcdf_gumbel(double p, double mu, double sigma) {
  if (sigma <= 0 || p < 0 || p > 1)
    return NAN;
  return mu - sigma * std::log(-std::log(p));
}


// [[Rcpp::export]]
NumericVector cpp_dgumbel(NumericVector x,
                          NumericVector mu, NumericVector sigma,
                          bool log_prob = false) {

  double z;
  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_gumbel(x[i % n], mu[i % nm], sigma[i % ns]);

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pgumbel(NumericVector x,
                          NumericVector mu, NumericVector sigma,
                          bool lower_tail = true, bool log_prob = false) {

  double z;
  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_gumbel(x[i % n], mu[i % nm], sigma[i % ns]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qgumbel(NumericVector p,
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
    q[i] = invcdf_gumbel(q[i % n], mu[i % nm], sigma[i % ns]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rgumbel(int n,
                          NumericVector mu, NumericVector sigma) {

  double u;
  int nm = mu.length();
  int ns = sigma.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    u = R::runif(0, 1);
    x[i] = invcdf_gumbel(u, mu[i % nm], sigma[i % ns]);
  }

  return x;
}

