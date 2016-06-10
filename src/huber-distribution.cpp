#include <Rcpp.h>
#include "shared.h"

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


double pdf_huber(double x, double mu, double sigma, double c) {
  if (sigma <= 0.0 || c <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  double z, A, rho;
  z = abs((x - mu)/sigma);
  A = 2.0*sqrt(2.0*M_PI) * (Phi(c) + phi(c)/c - 0.5);

  if (z <= c)
    rho = pow(z, 2.0)/2.0;
  else
    rho = c*z - pow(c, 2.0)/2.0;

  return exp(-rho)/A/sigma;
}

double cdf_huber(double x, double mu, double sigma, double c) {
  if (sigma <= 0.0 || c <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }

  double A, z, az, p;
  A = 2.0*(phi(c)/c - Phi(-c) + 0.5);
  z = (x - mu)/sigma;
  az = -abs(z);
  
  if (az <= -c) 
    p = exp(pow(c, 2.0)/2.0)/c * exp(c*az) / sqrt(2.0*M_PI)/A;
  else
    p = (phi(c)/c + Phi(az) - Phi(-c))/A;
  
  if (z <= 0.0)
    return p;
  else
    return 1.0 - p;
}

double invcdf_huber(double p, double mu, double sigma, double c) {
  if (sigma <= 0.0 || c <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  double x, pm, A;
  A = 2.0*sqrt(2.0*M_PI) * (Phi(c) + phi(c)/c - 0.5);
  pm = std::min(p, 1.0 - p);

  if (pm <= sqrt(2.0*M_PI) * phi(c)/(c*A))
    x = log(c*pm*A)/c - c/2.0;
  else
    x = InvPhi(abs(1.0 - Phi(c) + pm*A/sqrt(2.0*M_PI) - phi(c)/c));

  if (p < 0.5)
    return mu + x*sigma;
  else
    return mu - x*sigma;
}


// [[Rcpp::export]]
NumericVector cpp_dhuber(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& epsilon,
    bool log_prob = false
  ) {
  
  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int ne = epsilon.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, ne));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_huber(x[i % n], mu[i % nm], sigma[i % ns], epsilon[i % ne]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_phuber(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& epsilon,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int ne = epsilon.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, ne));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_huber(x[i % n], mu[i % nm], sigma[i % ns], epsilon[i % ne]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qhuber(
    const NumericVector& p,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& epsilon,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = p.length();
  int nm = mu.length();
  int ns = sigma.length();
  int ne = epsilon.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, ne));
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1.0 - pp[i];
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_huber(pp[i % n], mu[i % nm], sigma[i % ns], epsilon[i % ne]);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rhuber(
    const int n,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& epsilon
  ) {
  
  double u;
  int nm = mu.length();
  int ns = sigma.length();
  int ne = epsilon.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++) {
    u = R::runif(0.0, 1.0);
    x[i] = invcdf_huber(u, mu[i % nm], sigma[i % ns], epsilon[i % ne]);
  }
  
  return x;
}

