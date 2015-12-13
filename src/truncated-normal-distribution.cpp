#include <Rcpp.h>
#include "shared.h"
using namespace Rcpp;

/*
*  Truncated Normal distribution
*
*  Values:
*  x
*
*  Parameters:
*  mu
*  sigma > 0
*  a, b
*
*  z = (x-mu)/sigma
*
*  f(x)    = phi(z) / (Phi((b-mu)/sigma) - Phi((mu-a)/sigma))
*  F(x)    = (Phi(z) - Phi((mu-a)/sigma)) / (Phi((b-mu)/sigma) - Phi((a-mu)/sigma))
*  F^-1(p) = Phi^-1(Phi((mu-a)/sigma) + p * (Phi((b-mu)/sigma) - Phi((a-mu)/sigma)))
*
*  where phi() is PDF for N(0, 1) and Phi() is CDF for N(0, 1)
*
*/

const double pisq = sqrt(2*M_PI);

double pdf_tnorm(double x, double mu, double sigma, double a, double b) {
  if (sigma <= 0 || b < a)
    return NAN;
  double Phi_a, Phi_b;
  if (x > a && x < b) {
    Phi_a = Phi((a-mu)/sigma);
    Phi_b = Phi((b-mu)/sigma);
    return std::exp(-std::pow(x-mu, 2) / (2*std::pow(sigma, 2))) / (pisq*sigma * (Phi_b - Phi_a));
  } else {
    return 0;
  }
}

double cdf_tnorm(double x, double mu, double sigma, double a, double b) {
  if (sigma <= 0 || b < a)
    return NAN;
  double Phi_x, Phi_a, Phi_b;
  if (x > a && x < b) {
    Phi_x = Phi((x-mu)/sigma);
    Phi_a = Phi((a-mu)/sigma);
    Phi_b = Phi((b-mu)/sigma);
    return (Phi_x - Phi_a) / (Phi_b - Phi_a);
  } else if (x >= b) {
    return 1;
  } else {
    return 0;
  }
}

double invcdf_tnorm(double p, double mu, double sigma, double a, double b) {
  if (sigma <= 0 || b < a || p < 0 || p > 1)
    return NAN;
  double Phi_a, Phi_b;
  Phi_a = Phi((a-mu)/sigma);
  Phi_b = Phi((b-mu)/sigma);
  return InvPhi(Phi_a + p * (Phi_b - Phi_a));
}

double rng_tnorm(double mu, double sigma, double a, double b) {
  if (sigma <= 0 || b < a)
    return NAN;

  double r, u, za, zb;
  bool stop = false;

  za = (a-mu)/sigma;
  zb = (b-mu)/sigma;

  if (zb - za < pisq) {
    if (0 < za) {
      while (!stop) {
        r = R::runif(za, zb);
        u = R::runif(0, 1);
        stop = (u <= std::exp((std::pow(za, 2) - std::pow(r, 2))/2));
      }
    } else if (zb < 0) {
      while (!stop) {
        r = R::runif(za, zb);
        u = R::runif(0, 1);
        stop = (u <= std::exp((std::pow(zb, 2) - std::pow(r, 2))/2));
      }
    } else {
      while (!stop) {
        r = R::runif(za, zb);
        u = R::runif(0, 1);
        stop = (u <= std::exp(-std::pow(r, 2)/2));
      }
    }
  } else {
    while (!stop) {
      r = R::rnorm(0, 1);
      stop = (r > za && r < zb);
    }
  }

  return mu + sigma * r;
}


// [[Rcpp::export]]
NumericVector cpp_dtnorm(NumericVector x,
                         NumericVector mu, NumericVector sigma,
                         NumericVector a, NumericVector b,
                         bool log_prob = false) {

  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int na = mu.length();
  int nb = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, na, nb));
  NumericVector p(n);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_tnorm(x[i % n], mu[i % nm], sigma[i % ns], a[i % na], b[i % nb]);

  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = std::log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_ptnorm(NumericVector x,
                         NumericVector mu, NumericVector sigma,
                         NumericVector a, NumericVector b,
                         bool lower_tail = true, bool log_prob = false) {

  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int na = mu.length();
  int nb = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, na, nb));
  NumericVector p(n);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_tnorm(x[i % n], mu[i % nm], sigma[i % ns], a[i % na], b[i % nb]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qtnorm(NumericVector p,
                         NumericVector mu, NumericVector sigma,
                         NumericVector a, NumericVector b,
                         bool lower_tail = true, bool log_prob = false) {

  int n  = p.length();
  int nm = mu.length();
  int ns = sigma.length();
  int na = mu.length();
  int nb = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, na, nb));
  NumericVector q(n);

  for (int i = 0; i < Nmax; i++)
    q[i] = cdf_tnorm(p[i % n], mu[i % nm], sigma[i % ns], a[i % na], b[i % nb]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rtnorm(int n,
                         NumericVector mu, NumericVector sigma,
                         NumericVector a, NumericVector b) {

  int nm = mu.length();
  int ns = sigma.length();
  int na = mu.length();
  int nb = sigma.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_tnorm(mu[i % nm], sigma[i % ns], a[i % na], b[i % nb]);

  return x;
}

