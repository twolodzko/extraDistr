#include <Rcpp.h>
#include "shared.h"
using namespace Rcpp;

double pmf_bpois(double x, double y, double a, double b, double c) {
  
  if (a < 0 || b < 0 || c < 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (!isInteger(x))
    return 0;
  
  if (floor(y) != y) {
    char msg[55];
    int msg_len = sprintf(msg, "non-integer y = %f", y);
    if (msg_len >= 55 - 1 || msg_len < 0)
      Rcpp::warning("non-integer y");
    else
      Rcpp::warning(msg);
    Rcpp::warning(msg);
    return 0;
  }
  
  double tmp = exp(-(a+b+c)); 
  tmp *= (pow(a, x) / factorial(x)) * (pow(b, y) / factorial(y));
  double xy = 0;
  
  if (x < y) {
    for (int k = 0; k < x; k++)
      xy += R::choose(x, k) * R::choose(y, k) * factorial(k) * pow(c/(a*b), k);
  } else {
    for (int k = 0; k < y; k++)
      xy += R::choose(x, k) * R::choose(y, k) * factorial(k) * pow(c/(a*b), k);
  }
  
  return tmp * xy;
}


double logpmf_bpois(double x, double y, double a, double b, double c) {
  
  if (a < 0 || b < 0 || c < 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (!isInteger(x))
    return -INFINITY;
  
  if (floor(y) != y) {
    char msg[55];
    sprintf(msg, "non-integer y = %f", x);
    Rcpp::warning(msg);
    return -INFINITY;
  }
  
  double tmp = -(a+b+c); 
  tmp += log(pow(x, a) / lfactorial(x)) + log(pow(y, b) / lfactorial(y));
  double xy = 0;
  
  if (x < y) {
    for (int k = 0; k < x; k++)
      xy += R::choose(x, k) * R::choose(y, k) * lfactorial(k) * pow(c/(a*b), k);
  } else {
    for (int k = 0; k < y; k++)
      xy += R::choose(x, k) * R::choose(y, k) * lfactorial(k) * pow(c/(a*b), k);
  }
  
  return tmp + log(xy);
}



// [[Rcpp::export]]
NumericVector cpp_dbpois(
    const NumericVector& x,
    const NumericVector& y,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c,
    bool log_prob = false
  ) {
  
  int nx = x.length();
  int ny = y.length();
  int na = a.length();
  int nb = b.length();
  int nc = c.length();
  int Nmax = Rcpp::max(IntegerVector::create(nx, ny, na, nb, nc));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pmf_bpois(x[i % nx], y[i % ny], a[i % na], b[i % nb], c[i % nc]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}

