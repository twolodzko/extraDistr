#include <Rcpp.h>
using namespace Rcpp;


/*
 *  Pareto distribution
 *
 *  Values:
 *  x
 *
 *  Parameters:
 *  a, b > 0
 *
 *  f(x)    = (a*b^a) / x^{a+1}
 *  F(x)    = 1 - (b/x)^a
 *  F^-1(p) = b/(1-p)^{1-a}
 *
 */

double pdf_pareto(double x, double a, double b) {
  if (x >= b)
    return a * pow(b, a) / pow(x, a+1);
  else
    return 0;
}

double cdf_pareto(double x, double a, double b) {
  if (x >= b)
    return 1 - pow(b/x, a);
  else
    return 0;
}

double invcdf_pareto(double p, double a, double b) {
  return b / pow(1-p, 1/a);
}

double logpdf_pareto(double x, double a, double b) {
  if (x >= b)
    return log(a) + log(b)*a - log(x)*(a+1);
  else
    return -INFINITY;
}

double invcdf_pareto2(double p, double a, double b) {
  return exp(log(b) - log(1-p)*(1/a));
}


// [[Rcpp::export]]
NumericVector cpp_dpareto(NumericVector x,
                          NumericVector a, NumericVector b,
                          bool log_prob = false) {

  if (is_true(any(a <= 0)) || is_true(any(b <= 0)))
    throw Rcpp::exception("Values of a and b should be > 0.");

  int n = x.length();
  int na = a.length();
  int nb = b.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_pareto(x[i % n], a[i % na], b[i % nb]);

  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = exp(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_ppareto(NumericVector x,
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
    p[i] = cdf_pareto(x[i % n], a[i % na], b[i % nb]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qpareto(NumericVector p,
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
    q[i] = invcdf_pareto(p[i % n], a[i % na], b[i % nb]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rpareto(int n,
                          NumericVector a, NumericVector b) {

  if (is_true(any(a <= 0)) || is_true(any(b <= 0)))
    throw Rcpp::exception("Values of a and b should be > 0.");

  double u;
  int na = a.length();
  int nb = b.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    u = R::runif(0, 1);
    x[i] = invcdf_pareto(u, a[i % na], b[i % nb]);
  }

  return x;
}

