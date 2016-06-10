#include <Rcpp.h>

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


/*
*  Triangular distribution
*
*  Values:
*  x
*
*  Parameters:
*  a
*  b > a
*  a <= c <= b
*
*  f(x)    = { (2*(x-a)) / ((b-a)*(c-a))  x < c
*            { 2/(b-a)                    x = c
*            { (2*(b-x)) / ((b-a)*(b-c))  x > c
*  F(x)    = { (x-a)^2 / ((b-a)*(c-a))
*            { 1 - ((b-x)^2 / ((b-a)*(b-c)))
*  F^-1(p) = { a + sqrt(p*(b-a)*(c-a))    p < (c-a)/(b-a)
*            { b - sqrt((1-p)*(b-a)*(b-c));
*/

double pdf_triangular(double x, double a, double b, double c) {
  if (a > c || c > b || a == b) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < a || x > b) {
    return 0.0;
  } else if (x < c) {
    return 2.0*(x-a) / ((b-a)*(c-a));
  } else if (x > c) {
    return 2.0*(b-x) / ((b-a)*(b-c));
  } else {
    return 2.0/(b-a);
  }
}

double cdf_triangular(double x, double a, double b, double c) {
  if (a > c || c > b || a == b) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < a) {
    return 0.0;
  } else if (x >= b) {
    return 1;
  } else if (x <= c) {
    return pow(x-a, 2.0) / ((b-a)*(c-a));
  } else {
    return 1.0 - (pow(b-x, 2.0) / ((b-a)*(b-c)));
  }
}

double invcdf_triangular(double p, double a, double b, double c) {
  if (a > c || c > b || a == b || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double fc = (c-a)/(b-a);
  if (p < fc) {
    return a + sqrt(p*(b-a)*(c-a));
  } else {
    return b - sqrt((1.0-p)*(b-a)*(b-c));
  }
}

double rng_triangular(double a, double b, double c) {
  if (a > c || c > b || a == b) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double u, v, r, cc;
  r = b - a;
  cc = (c-a)/r;
  u = R::runif(0.0, 1.0);
  v = R::runif(0.0, 1.0);
  return ((1.0-cc) * std::min(u, v) + cc * std::max(u, v)) * r + a;
}


// [[Rcpp::export]]
NumericVector cpp_dtriang(
    const NumericVector& x,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c,
    bool log_prob = false
  ) {

  int n = x.length();
  int na = a.length();
  int nb = b.length();
  int nc = c.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb, nc));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_triangular(x[i % n], a[i % na], b[i % nb], c[i % nc]);

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_ptriang(
    const NumericVector& x,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c,
    bool lower_tail = true, bool log_prob = false
  ) {

  int n  = x.length();
  int na = a.length();
  int nb = b.length();
  int nc = c.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb, nc));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_triangular(x[i % n], a[i % na], b[i % nb], c[i % nc]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qtriang(
    const NumericVector& p,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c,
    bool lower_tail = true, bool log_prob = false
  ) {

  int n  = p.length();
  int na = a.length();
  int nb = b.length();
  int nc = c.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb, nc));
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);

  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1.0 - pp[i];

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_triangular(pp[i % n], a[i % na], b[i % nb], c[i % nc]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rtriang(
    const int n,
    const NumericVector& a,
    const NumericVector& b,
    const NumericVector& c
  ) {

  int na = a.length();
  int nb = b.length();
  int nc = c.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_triangular(a[i % na], b[i % nb], c[i % nc]);

  return x;
}

