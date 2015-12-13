#include <Rcpp.h>
#include "shared.h"
using namespace Rcpp;

double pmf_bpois(double x, double y, double a, double b, double c) {
  
  if (std::floor(x) != x || std::floor(y) != y)
    return 0;
  
  if (a < 0 || b < 0 || c < 0)
    return NAN;
  
  double tmp = std::exp(-(a+b+c)); 
  tmp *= (std::pow(a, x) / factorial(x)) * (std::pow(b, y) / factorial(y));
  double xy = 0;
  
  if (x < y) {
    for (int k = 0; k < x; k++)
      xy += R::choose(x, k) * R::choose(y, k) * factorial(k) * std::pow(c/(a*b), k);
  } else {
    for (int k = 0; k < y; k++)
      xy += R::choose(x, k) * R::choose(y, k) * factorial(k) * std::pow(c/(a*b), k);
  }
  
  return tmp * xy;
}


double logpmf_bpois(double x, double y, double a, double b, double c) {
  
  if (std::floor(x) != x || std::floor(y) != y)
    return -INFINITY;
  
  if (a < 0 || b < 0 || c < 0)
    return NAN;
  
  double tmp = -(a+b+c); 
  tmp += std::log(std::pow(x, a) / lfactorial(x)) + std::log(pow(y, b) / lfactorial(y));
  double xy = 0;
  
  if (x < y) {
    for (int k = 0; k < x; k++)
      xy += R::choose(x, k) * R::choose(y, k) * lfactorial(k) * std::pow(c/(a*b), k);
  } else {
    for (int k = 0; k < y; k++)
      xy += R::choose(x, k) * R::choose(y, k) * lfactorial(k) * std::pow(c/(a*b), k);
  }
  
  return tmp + std::log(xy);
}



// [[Rcpp::export]]
NumericVector cpp_dbpois(NumericVector x, NumericVector y,
                         NumericVector a, NumericVector b, NumericVector c,
                         bool log_prob = false) {
  
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
      p[i] = std::log(p[i]);
  
  return p;
}

