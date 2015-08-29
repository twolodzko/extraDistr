#include <Rcpp.h>
using namespace Rcpp;


/*
*  Gompertz distribution
*
*  Values:
*  x >= 0
*
*  Parameters:
*  a > 0
*  b > 0
*
*  f(x)    = a*exp(b*x - a/b * (exp(bx)-1))
*  F(x)    = 1-exp(-a/b * (exp(bx)-1))
*  F^-1(p) = 1/b * log(1 - b/a * log(1-p))
*
* References:
*
* Lenart, A. (2012). The Gompertz distribution and Maximum Likelihood Estimation
* of its parameters - a revision. MPIDR WORKING PAPER WP 2012-008.
*
*/


double pdf_gompertz(double x, double a, double b) {
  if (x >= 0)
    return a * exp(b*x - a/b * (exp(b*x) - 1));
  else
    return 0;
}

double cdf_gompertz(double x, double a, double b) {
  if (x >= 0)
    return 1 - exp(-a/b * (exp(b*x) - 1));
  else
    return 0;
}

double invcdf_gompertz(double p, double a, double b) {
  return 1/b * log(1 - b/a * log(1-p));
}

double logpdf_gompertz(double x, double a, double b) {
  if (x >= 0)
    return log(a) + (b*x - a/b * (exp(b*x) - 1));
  else
    return -INFINITY;
}


// [[Rcpp::export]]
NumericVector cpp_dgompertz(NumericVector x,
                            NumericVector a, NumericVector b,
                            bool log_prob = false) {

  if (is_true(any(a <= 0)) || is_true(any(b <= 0)))
    throw Rcpp::exception("Values of a and b should be > 0.");

  int n  = x.length();
  int na = a.length();
  int nb = b.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_gompertz(x[i % n], a[i % na], b[i % nb]);

  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = exp(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pgompertz(NumericVector x,
                            NumericVector a, NumericVector b,
                            bool lower_tail = true, bool log_prob = false) {

  if (is_true(any(a <= 0)) || is_true(any(b <= 0)))
    throw Rcpp::exception("Values of a and b should be > 0.");

  int n  = x.length();
  int na = a.length();
  int nb = b.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_gompertz(x[i % n], a[i % na], b[i % nb]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qgompertz(NumericVector p,
                            NumericVector a, NumericVector b,
                            bool lower_tail = true, bool log_prob = false) {

  if (is_true(any(p < 0)) || is_true(any(p > 1)))
    throw Rcpp::exception("Probabilities should range from 0 to 1.");
  if (is_true(any(a <= 0)) || is_true(any(b <= 0)))
    throw Rcpp::exception("Values of a and b should be > 0.");

  int n  = p.length();
  int na = a.length();
  int nb = b.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector q(Nmax);

  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = exp(p[i]);

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      p[i] = 1-p[i];

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_gompertz(p[i % n], a[i % na], b[i % nb]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rgompertz(int n,
                            NumericVector a, NumericVector b) {

  if (is_true(any(a <= 0)) || is_true(any(b <= 0)))
    throw Rcpp::exception("Values of a and b should be > 0.");

  double u;
  int na = a.length();
  int nb = b.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    u = R::runif(0, 1);
    x[i] = invcdf_gompertz(u, a[i % na], b[i % nb]);
  }

  return x;
}

