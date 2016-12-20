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
  
  if (!isInteger(x) || x < 0.0 || !R_finite(x) || !R_finite(y))
    return 0.0;
  
  if (!isInteger(y, false)) {
    char msg[55];
    std::snprintf(msg, sizeof(msg), "non-integer y = %f", y);
    Rcpp::warning(msg);
    return 0.0;
  }
  
  if (y < 0.0)
    return 0.0;
  
  double tmp = exp(-(a+b+c)); 
  tmp *= (pow(a, x) / factorial(x)) * (pow(b, y) / factorial(y));
  
  double z = (x < y) ? x : y;
  double c_ab = c/(a*b);
  double xy = 0.0;
  
  for (double k = 0.0; k <= z; k += 1.0)
    xy += R::choose(x, k) * R::choose(y, k) * factorial(k) * pow(c_ab, k);
  
  return tmp * xy;
}


// [[Rcpp::export]]
NumericVector cpp_dbpois(
    const NumericVector& x,
    const NumericVector& y,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(y.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  dims.push_back(c.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  if (dims[0] != dims[1])
    Rcpp::stop("lengths of x and y differ");
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pmf_bpois(x[i % dims[0]], y[i % dims[1]],
                     a[i % dims[2]], b[i % dims[3]], c[i % dims[4]]);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericMatrix cpp_rbpois(
    const int& n,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c
  ) {
  
  double u, v, w;
  std::vector<int> dims;
  dims.push_back(a.length());
  dims.push_back(b.length());
  dims.push_back(c.length());
  NumericMatrix x(n, 2);
  
  for (int i = 0; i < n; i++) {
    if (ISNAN(a[i % dims[0]]) || ISNAN(b[i % dims[1]]) || ISNAN(c[i % dims[2]])) {
      x(i, 0) = NA_REAL;
      x(i, 1) = NA_REAL;
    } else if (a[i % dims[0]] < 0.0 || b[i % dims[1]] < 0.0 || c[i % dims[2]] < 0.0) {
      Rcpp::warning("NaNs produced");
      x(i, 0) = NAN;
      x(i, 1) = NAN;
    } else {
      u = R::rpois(a[i % dims[0]]);
      v = R::rpois(b[i % dims[1]]);
      w = R::rpois(c[i % dims[2]]);
      x(i, 0) = u+w;
      x(i, 1) = v+w;
    }
  }

  return x;
}

