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


double pmf_bpois(double x, double y, double a, double b, double c) {
  
  if (ISNAN(x) || ISNAN(y) || ISNAN(a) || ISNAN(b) || ISNAN(c))
    return NA_REAL;
  
  if (a < 0.0 || b < 0.0 || c < 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (!isInteger(x) || x < 0.0)
    return 0.0;
  
  if (floor(y) != y) {
    char msg[55];
    std::snprintf(msg, sizeof(msg), "non-integer y = %f", y);
    Rcpp::warning(msg);
    return 0.0;
  }
  
  if (y < 0.0)
    return 0.0;
  
  double tmp = exp(-(a+b+c)); 
  tmp *= (pow(a, x) / factorial(x)) * (pow(b, y) / factorial(y));
  double xy = 0.0;
  
  double z;
  if (x < y)
    z = x;
  else
    z = y;
  
  double k = 0.0;
  while (k < z) {
    xy += R::choose(x, k) * R::choose(y, k) * factorial(k) * pow(c/(a*b), k);
    k += 1.0;
  }
  
  return tmp * xy;
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
    if (ISNAN(a[i % na]) || ISNAN(b[i % nb]) || ISNAN(c[i % nc])) {
      x(i, 0) = NA_REAL;
      x(i, 1) = NA_REAL;
    } else if (a[i % na] < 0.0 || b[i % nb] < 0.0 || c[i % nc] < 0.0) {
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

