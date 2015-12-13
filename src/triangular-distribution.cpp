#include <Rcpp.h>
using namespace Rcpp;


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
  if (a > c || c > b)
    return NAN;
  if (x < a || x > b) {
    return 0;
  } else if (x < c) {
    return 2*(x-a) / ((b-a)*(c-a));
  } else if (x > c) {
    return 2*(b-x) / ((b-a)*(b-c));
  } else {
    return 2/(b-a);
  }
}

double cdf_triangular(double x, double a, double b, double c) {
  if (a > c || c > b)
    return NAN;
  if (x < a) {
    return 0;
  } else if (x >= b) {
    return 1;
  } else if (x <= c) {
    return std::pow(x-a, 2) / ((b-a)*(c-a));
  } else {
    return 1 - (std::pow(b-x, 2) / ((b-a)*(b-c)));
  }
}

double invcdf_triangular(double p, double a, double b, double c) {
  if (a > c || c > b || p < 0 || p > 1)
    return NAN;
  double fc = (c-a)/(b-a);
  if (p < fc) {
    return a + std::sqrt(p*(b-a)*(c-a));
  } else {
    return b - std::sqrt((1-p)*(b-a)*(b-c));
  }
}


// [[Rcpp::export]]
NumericVector cpp_dtriang(NumericVector x,
                          NumericVector a, NumericVector b, NumericVector c,
                          bool log_prob = false) {

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
      p[i] = std::log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_ptriang(NumericVector x,
                          NumericVector a, NumericVector b, NumericVector c,
                          bool lower_tail = true, bool log_prob = false) {

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
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qtriang(NumericVector p,
                          NumericVector a, NumericVector b, NumericVector c,
                          bool lower_tail = true, bool log_prob = false) {

  int n  = p.length();
  int na = a.length();
  int nb = b.length();
  int nc = c.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb, nc));
  NumericVector q(Nmax);

  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = std::exp(p[i]);

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      p[i] = 1-p[i];

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_triangular(p[i % n], a[i % na], b[i % nb], c[i % nc]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rtriang(int n,
                          NumericVector a, NumericVector b, NumericVector c) {

  double u;
  int na = a.length();
  int nb = b.length();
  int nc = c.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    u = R::runif(0, 1);
    x[i] = invcdf_triangular(u, a[i % na], b[i % nb], c[i % nc]);
  }

  return x;
}

