#include <Rcpp.h>
#include "shared.h"
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;


inline double logpmf_bpois(double x, double y, double a, double b, double c,
                           bool& throw_warning) {
  
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(y) || ISNAN(a) || ISNAN(b) || ISNAN(c))
    return x+y+a+b+c;
#endif
  
  if (a < 0.0 || b < 0.0 || c < 0.0) {
    throw_warning = true;
    return NAN;
  }
  
  if (!isInteger(x) || x < 0.0 || !R_FINITE(x) ||
      !R_FINITE(y) || !isInteger(y)) {
      return R_NegInf;
  }
  
  if (y < 0.0)
    return R_NegInf;
  
  // exp(-(a+b+c))
  double tmp = -(a+b+c); 
  // tmp *= (pow(a, x) / factorial(x)) * (pow(b, y) / factorial(y));
  tmp += (log(a) * x - lfactorial(x)) + (log(b) * y - lfactorial(y));
  
  double minxy = static_cast<int>( (x < y) ? x : y );
  // c_ab = c/(a*b)
  double lc_ab = log(c) - log(a) - log(b);
  
  double dk;
  double mx = R_NegInf;
  std::vector<double> ls(minxy+1);
  
  for (int k = 0; k <= minxy; k++) {
    dk = static_cast<double>(k);
    // xy += R::choose(x, k) * R::choose(y, k) * factorial(k) * pow(c_ab, k);
    ls[k] = R::lchoose(x, dk) + R::lchoose(y, dk) + lfactorial(dk) + lc_ab * dk;
    if (ls[k] > mx)
      mx = ls[k];
  }
  
  double xy = 0.0;
  
  for (int k = 0; k <= minxy; k++)
    xy += exp(ls[k] - mx);    // log-sum-exp trick
  
  xy = log(xy) + mx;
  
  return tmp + xy;
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
  
  if (std::min({x.length(), y.length(),
                a.length(), b.length(),
                c.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    y.length(),
    a.length(),
    b.length(),
    c.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  if (x.length() != y.length())
    Rcpp::stop("lengths of x and y differ");
  
  for (int i = 0; i < Nmax; i++)
    p[i] = logpmf_bpois(GETV(x, i), GETV(y, i), GETV(a, i),
                        GETV(b, i), GETV(c, i), throw_warning);
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
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
  
  if (std::min({a.length(), b.length(), c.length()}) < 1) {
    Rcpp::warning("NAs produced");
    NumericMatrix out(n, 2);
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }
  
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

