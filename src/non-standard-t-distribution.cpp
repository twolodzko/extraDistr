#include <Rcpp.h>
using namespace Rcpp;


/*
*  Non-standard t-distribution
*
*  Values:
*  x
*
*  Parameters:
*  df > 0
*  mu
*  sigma > 0
*
*/

double pdf_nst(double x, double df, double mu, double sigma) {
  if (df <= 0 || sigma <= 0)
    return NAN;
  if (df == 1)
    return R::dcauchy(x, mu, sigma, false);
  double z = (x - mu)/sigma;
  return R::dt(z, df, false)/sigma;
}

double cdf_nst(double x, double df, double mu, double sigma) {
  if (df <= 0 || sigma <= 0)
    return NAN;
  if (df == 1)
    return R::pcauchy(x, mu, sigma, true, false);
  double z = (x - mu)/sigma;
  return R::pt(z, df, true, false);
}

double invcdf_nst(double p, double df, double mu, double sigma) {
  if (df <= 0 || sigma <= 0)
    return NAN;
  if (df == 1)
    return R::qcauchy(p, mu, sigma, true, false);
  return R::qt(p, df, true, false)*sigma + mu;
}

double rng_nst(double df, double mu, double sigma) {
  if (df <= 0 || sigma <= 0)
    return NAN;
  if (df == 1)
    return R::rcauchy(mu, sigma);
  return R::rt(df)*sigma + mu;
}


// [[Rcpp::export]]
NumericVector cpp_dnst(NumericVector x,
                       NumericVector df, NumericVector mu, NumericVector sigma,
                       bool log_prob = false) {
  
  double z;
  int n  = x.length();
  int nd = df.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nd, nm, ns));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_nst(x[i % n], df[i % nd], mu[i % nm], sigma[i % ns]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pnst(NumericVector x,
                       NumericVector df, NumericVector mu, NumericVector sigma,
                       bool lower_tail = true, bool log_prob = false) {
  
  double z;
  int n  = x.length();
  int nd = df.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nd, nm, ns));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_nst(x[i % n], df[i % nd], mu[i % nm], sigma[i % ns]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qnst(NumericVector p,
                       NumericVector df, NumericVector mu, NumericVector sigma,
                       bool lower_tail = true, bool log_prob = false) {
  
  int n  = p.length();
  int nd = df.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nd, nm, ns));
  NumericVector q(Nmax);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = std::exp(p[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      p[i] = 1-p[i];
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_nst(p[i % n], df[i % nd], mu[i % nm], sigma[i % ns]);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rnst(int n,
                       NumericVector df, NumericVector mu, NumericVector sigma) {
  
  double u;
  int nd = df.length();
  int nm = mu.length();
  int ns = sigma.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_nst(df[i % nd], mu[i % nm], sigma[i % ns]);
  
  return x;
}

