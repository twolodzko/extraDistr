#include <Rcpp.h>

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using std::sin;
using std::cos;
using std::tan;
using std::atan;
using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;


/*
*  Non-standard t-distribution
*
*  Values:
*  x
*
*  Parameters:
*  nu > 0
*  mu
*  sigma > 0
*
*/

double pdf_nst(double x, double nu, double mu, double sigma) {
  if (nu <= 0 || sigma <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (nu == 1)
    return R::dcauchy(x, mu, sigma, false);
  double z = (x - mu)/sigma;
  return R::dt(z, nu, false)/sigma;
}

double cdf_nst(double x, double nu, double mu, double sigma) {
  if (nu <= 0 || sigma <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (nu == 1)
    return R::pcauchy(x, mu, sigma, true, false);
  double z = (x - mu)/sigma;
  return R::pt(z, nu, true, false);
}

double invcdf_nst(double p, double nu, double mu, double sigma) {
  if (nu <= 0 || sigma <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (nu == 1)
    return R::qcauchy(p, mu, sigma, true, false);
  return R::qt(p, nu, true, false)*sigma + mu;
}

double rng_nst(double nu, double mu, double sigma) {
  if (nu <= 0 || sigma <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (nu == 1)
    return R::rcauchy(mu, sigma);
  return R::rt(nu)*sigma + mu;
}


// [[Rcpp::export]]
NumericVector cpp_dnst(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma,
    bool log_prob = false
  ) {
  
  int n  = x.length();
  int nn = nu.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, nm, ns));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_nst(x[i % n], nu[i % nn], mu[i % nm], sigma[i % ns]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pnst(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = x.length();
  int nn = nu.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, nm, ns));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_nst(x[i % n], nu[i % nn], mu[i % nm], sigma[i % ns]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qnst(
    const NumericVector& p,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = p.length();
  int nn = nu.length();
  int nm = mu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, nm, ns));
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1-pp[i];
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_nst(pp[i % n], nu[i % nn], mu[i % nm], sigma[i % ns]);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rnst(
    const int n,
    const NumericVector& nu,
    const NumericVector& mu,
    const NumericVector& sigma
  ) {
  
  int nn = nu.length();
  int nm = mu.length();
  int ns = sigma.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_nst(nu[i % nn], mu[i % nm], sigma[i % ns]);
  
  return x;
}

