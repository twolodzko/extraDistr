#include <Rcpp.h>
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


// [[Rcpp::export]]
NumericVector cpp_dlgser(IntegerVector x,
                         NumericVector theta,
                         bool log_prob = false) {

  if (is_true(any(theta <= 0)) || is_true(any(theta >= 1)))
    throw Rcpp::exception("Values of theta should fit 0 < theta < 1.");

  int n = x.length();
  int nt = theta.length();
  int Nmax = std::max(n, nt);
  NumericVector p(Nmax), a(nt);

  for (int j = 0; j < nt; j++)
    a[j] = -1/log(1-theta[j]);

  for (int i = 0; i < Nmax; i++) {
    if (x[i % n] >= 1)
      p[i] = a[i % nt] * pow(theta[i % nt], x[i % n]) / x[i % n];
    else
      p[i] = 0;
  }

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_plgser(IntegerVector x,
                         NumericVector theta,
                         bool lower_tail = true, bool log_prob = false) {

  if (is_true(any(theta <= 0)) || is_true(any(theta >= 1)))
    throw Rcpp::exception("Values of theta should fit 0 < theta < 1.");

  int n = x.length();
  int nt = theta.length();
  int Nmax = std::max(n, nt);
  NumericVector p(Nmax), a(nt);
  double b;

  for (int j = 0; j < nt; j++)
    a[j] = -1/log(1-theta[j]);

  for (int i = 0; i < Nmax; i++) {
    if (x[i % n] >= 1) {
      b = 0;
      for (int k = 1; k < x[i % n]+1; k++)
        b += pow(theta[i % nt], k) / k;
      p[i] = a[i % nt] * b;
    } else {
      p[i] = 0;
    }
  }

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rlgser(int n,
                         NumericVector theta) {

  if (is_true(any(theta <= 0)) || is_true(any(theta >= 1)))
    throw Rcpp::exception("Values of theta should fit 0 < theta < 1.");

  int nt = theta.length();
  NumericVector x(n), a(nt);
  double pk, k, u;

  for (int j = 0; j < nt; j++)
    a[j] = -theta[j  % nt] / log(1-theta[j]);

  pk = a[0];

  for (int i = 0; i < n; i++) {
    u = R::runif(0, 1);
    k = 1;
    while (u > pk) {
      u -= pk;
      pk = pk * theta[i % nt] * k/(k+1);
      k++;
    }
    x[i] = k;
    pk = a[i % nt];
  }

  return x;
}

