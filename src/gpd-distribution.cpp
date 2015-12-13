#include <Rcpp.h>
using namespace Rcpp;


/*
*  Generalized Pareto distribution
*
*  Values:
*  x
*
*  Parameters:
*  mu
*  sigma > 0
*  xi
*
*  z = (x-mu)/sigma
*  where 1+xi*z > 0
*
*  f(x)    = { (1+xi*z)^{-(xi+1)/xi}             if xi != 0
*            { exp(-z)                           otherwise
*  F(x)    = { 1-(1+xi*z)^{-1/xi}                if xi != 0
*            { 1-exp(-z)                         otherwise
*  F^-1(p) = { mu + sigma * ((1-p)^{-xi}-1)/xi   if xi != 0
*            { mu - sigma * log(1-p)             otherwise
*
*/

double pdf_gpd(double x, double mu, double sigma, double xi) {
  if (sigma <= 0)
    return NAN;
  double z = (x-mu)/sigma;
  if (xi != 0) {
    if (x >= mu)
      return 1-std::pow(1+xi*z, -1/xi);
    else
      return 0;
  } else {
    if (x >= mu && x <= mu - sigma/xi)
      return std::exp(-z);
    else
      return 0;
  }
}

double cdf_gpd(double x, double mu, double sigma, double xi) {
  if (sigma <= 0)
    return NAN;
  double z = (x-mu)/sigma;
  if (xi != 0) {
    if (x >= mu)
      return 1-std::pow(1+xi*z, -1/xi);
    else
      return 0;
  } else {
    if (x >= mu && x <= mu - sigma/xi)
      return 1-std::exp(-z);
    else
      return 0;
  }
}

double invcdf_gpd(double p, double mu, double sigma, double xi) {
  if (sigma <= 0 || p < 0 || p > 1)
    return NAN;
  if (xi != 0)
    return mu + sigma * (std::pow(1-p, -xi)-1)/xi;
  else
    return mu - sigma * std::log(1-p);
}


// [[Rcpp::export]]
NumericVector cpp_dgpd(NumericVector x,
                       NumericVector mu, NumericVector sigma, NumericVector xi,
                       bool log_prob = false) {

  double z;
  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int nx = xi.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, nx));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_gpd(x[i % n], mu[i % nm], sigma[i % ns], xi[i % nx]);

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pgpd(NumericVector x,
                       NumericVector mu, NumericVector sigma, NumericVector xi,
                       bool lower_tail = true, bool log_prob = false) {

  double z;
  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int nx = xi.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, nx));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_gpd(x[i % n], mu[i % nm], sigma[i % ns], xi[i % nx]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qgpd(NumericVector p,
                       NumericVector mu, NumericVector sigma, NumericVector xi,
                       bool lower_tail = true, bool log_prob = false) {

  int n  = p.length();
  int nm = mu.length();
  int ns = sigma.length();
  int nx = xi.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, nx));
  NumericVector q(Nmax);

  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = std::exp(p[i]);

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      p[i] = 1-p[i];

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_gpd(p[i % n], mu[i % nm], sigma[i % ns], xi[i % nx]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rgpd(int n,
                       NumericVector mu, NumericVector sigma, NumericVector xi) {

  double u;
  int nm = mu.length();
  int ns = sigma.length();
  int nx = xi.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    u = R::runif(0, 1);
    x[i] = invcdf_gpd(u, mu[i % nm], sigma[i % ns], xi[i % nx]);
  }

  return x;
}

