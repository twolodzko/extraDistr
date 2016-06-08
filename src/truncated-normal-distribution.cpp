#include <Rcpp.h>
#include "namespace.h"
#include "shared.h"


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

const double SQRTPI = sqrt(2*M_PI);

double pdf_tnorm(double x, double mu, double sigma, double a, double b) {
  if (sigma <= 0 || b <= a) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (a == -INFINITY && b == INFINITY)
    return R::dnorm(x, mu, sigma, false);
  
  double Phi_a, Phi_b;
  if (x > a && x < b) {
    Phi_a = Phi((a-mu)/sigma);
    Phi_b = Phi((b-mu)/sigma);
    return exp(-pow(x-mu, 2) / (2*pow(sigma, 2))) / (SQRTPI*sigma * (Phi_b - Phi_a));
  } else {
    return 0;
  }
}

double cdf_tnorm(double x, double mu, double sigma, double a, double b) {
  if (sigma <= 0 || b <= a) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (a == -INFINITY && b == INFINITY)
    return R::pnorm(x, mu, sigma, true, false);
  
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
  if (sigma <= 0 || b <= a || p < 0 || p > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (a == -INFINITY && b == INFINITY)
    return R::qnorm(p, mu, sigma, true, false);
  
  double Phi_a, Phi_b;
  Phi_a = Phi((a-mu)/sigma);
  Phi_b = Phi((b-mu)/sigma);
  return InvPhi(Phi_a + p * (Phi_b - Phi_a));
}

double rng_tnorm(double mu, double sigma, double a, double b) {
  if (sigma <= 0 || b <= a) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (a == -INFINITY && b == INFINITY)
    return R::rnorm(mu, sigma);

  double r, u, za, zb;
  bool stop = false;

  za = (a-mu)/sigma;
  zb = (b-mu)/sigma;

  if (zb - za < SQRTPI) {
    if (0 < za) {
      while (!stop) {
        r = R::runif(za, zb);
        u = R::runif(0, 1);
        stop = (u <= exp((pow(za, 2) - pow(r, 2))/2));
      }
    } else if (zb < 0) {
      while (!stop) {
        r = R::runif(za, zb);
        u = R::runif(0, 1);
        stop = (u <= exp((pow(zb, 2) - pow(r, 2))/2));
      }
    } else {
      while (!stop) {
        r = R::runif(za, zb);
        u = R::runif(0, 1);
        stop = (u <= exp(-pow(r, 2)/2));
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
NumericVector cpp_dtnorm(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& a,
    const NumericVector& b,
    bool log_prob = false
  ) {

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
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_ptnorm(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& a,
    const NumericVector& b,
    bool lower_tail = true, bool log_prob = false
  ) {

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
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qtnorm(
    const NumericVector& p,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& a,
    const NumericVector& b,
    bool lower_tail = true, bool log_prob = false
  ) {

  int n  = p.length();
  int nm = mu.length();
  int ns = sigma.length();
  int na = mu.length();
  int nb = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, na, nb));
  NumericVector q(n);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1-pp[i];

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_tnorm(pp[i % n], mu[i % nm], sigma[i % ns], a[i % na], b[i % nb]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rtnorm(
    const int n,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& a,
    const NumericVector& b
  ) {

  int nm = mu.length();
  int ns = sigma.length();
  int na = mu.length();
  int nb = sigma.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_tnorm(mu[i % nm], sigma[i % ns], a[i % na], b[i % nb]);

  return x;
}

