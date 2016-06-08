#include <Rcpp.h>
#include "namespace.h"
#include "shared.h"


double pmf_bpois(double x, double y, double a, double b, double c) {
  
  if (a < 0 || b < 0 || c < 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (!isInteger(x))
    return 0;
  
  if (floor(y) != y) {
    char msg[55];
    std::snprintf(msg, sizeof(msg), "non-integer y = %f", y);
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
    std::snprintf(msg, sizeof(msg), "non-integer y = %f", y);
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


// [[Rcpp::export]]
NumericMatrix cpp_rbpois(
    const int n,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c
  ) {
  
  double u, v, w;
  int na = a.length();
  int nb = b.length();
  int nc = c.length();
  NumericMatrix x(n, 2);
  
  for (int i = 0; i < n; i++) {
    if (a[i % na] < 0 || b[i % nb] < 0 || c[i % nc] < 0) {
      Rcpp::warning("NaNs produced");
      x(i, 0) = NAN;
      x(i, 1) = NAN;
    } else {
      u = R::rpois(a[i % na]);
      v = R::rpois(b[i % nb]);
      w = R::rpois(c[i % nc]);
      x(i, 0) = u+w;
      x(i, 1) = v+w;
    }
  }

  return x;
}

