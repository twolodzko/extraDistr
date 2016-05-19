#include <Rcpp.h>
#include "shared.h"
using namespace Rcpp;


/*
*  Logarithmic Series distribution
*
*  Values:
*  x
*
*  Parameters:
*  0 < theta < 1
*
*  f(x) = (-1/log(1-theta)*theta^x) / x
*  F(x) = -1/log(1-theta) * sum((theta^x)/x)
*
*/


double pdf_lgser(double x, double theta) {
  if (theta <= 0 || theta >= 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(x) || x < 1)
    return 0;
  double a = -1/log(1-theta);
  return a * pow(theta, x) / x;
}


double cdf_lgser(double x, double theta) {
  if (theta <= 0 || theta >= 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 1)
    return 0;
  if (std::isinf(x))
    return 1;
  
  double a = -1/log(1-theta);
  double b = 0;
  
  for (int k = 1; k < x+1; k++)
    b += pow(theta, k) / k;
  
  return a * b;
}


double invcdf_lgser(double p, double theta) {
  if (theta <= 0 || theta >= 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (p == 1)
    return INFINITY;
  
  double pk = -theta/log(1-theta);
  int k = 1;
  
  while (p > pk) {
    p -= pk;
    pk *= theta * k/(k+1);
    k++;
  }
  return k;
}


// [[Rcpp::export]]
NumericVector cpp_dlgser(
    const NumericVector& x,
    const NumericVector& theta,
    bool log_prob = false
  ) {

  int n = x.length();
  int nt = theta.length();
  int Nmax = std::max(n, nt);
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_lgser(x[i % n], theta[i % nt]);
 
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_plgser(
    const NumericVector& x,
    const NumericVector& theta,
    bool lower_tail = true, bool log_prob = false
  ) {

  int n = x.length();
  int nt = theta.length();
  int Nmax = std::max(n, nt);
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_lgser(x[i % n], theta[i % nt]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qlgser(
    const NumericVector& p,
    const NumericVector& theta,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n = p.length();
  int nt = theta.length();
  int Nmax = std::max(n, nt);
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1-pp[i];
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_lgser(pp[i % n], theta[i % nt]);
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rlgser(
    const int n,
    const NumericVector& theta
  ) {

  int nt = theta.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    double u = R::runif(0, 1);
    x[i] = invcdf_lgser(u, theta[i % nt]);
  }

  return x;
}

