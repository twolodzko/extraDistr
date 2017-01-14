#include <Rcpp.h>
#include "shared.h"

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;


inline double pmf_bpois(double x, double y, double a, double b, double c,
                        bool& throw_warning) {
  
  if (ISNAN(x) || ISNAN(y) || ISNAN(a) || ISNAN(b) || ISNAN(c))
    return x+y+a+b+c;

  if (a < 0.0 || b < 0.0 || c < 0.0) {
    throw_warning = true;
    return NAN;
  }
  
  if (!isInteger(x) || x < 0.0 || !R_FINITE(x) || !R_FINITE(y))
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
  
  bool throw_warning = false;
  
  if (dims[0] != dims[1])
    Rcpp::stop("lengths of x and y differ");
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pmf_bpois(GETV(x, i), GETV(y, i), GETV(a, i),
                     GETV(b, i), GETV(c, i), throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericMatrix cpp_rbpois(
    const int& n,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c
  ) {
  
  NumericMatrix x(n, 2);
  double u, v, w;
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++) {
    if (ISNAN(GETV(a, i)) || ISNAN(GETV(b, i)) || ISNAN(GETV(c, i)) || 
        GETV(a, i) < 0.0 || GETV(b, i) < 0.0 || GETV(c, i) < 0.0) {
      throw_warning = true;
      x(i, 0) = NA_REAL;
      x(i, 1) = NA_REAL;
    } else {
      u = R::rpois(GETV(a, i));
      v = R::rpois(GETV(b, i));
      w = R::rpois(GETV(c, i));
      x(i, 0) = u+w;
      x(i, 1) = v+w;
    }
  }

  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

