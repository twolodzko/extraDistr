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
 *  Rayleigh distribution
 *
 *  Values:
 *  x >= 0
 *
 *  Parameters:
 *  sigma > 0
 *
 *  f(x)    = x/sigma^2 * exp(-(x^2 / 2*sigma^2))
 *  F(x)    = 1 - exp(-x^2 / 2*sigma^2)
 *  F^-1(p) = sigma * sqrt(-2 * log(1-p))
 *
 */

double pdf_rayleigh(double x, double sigma) {
  if (sigma <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0 || std::isinf(x))
    return 0;
  return x/pow(sigma, 2.0) * exp(-pow(x, 2.0) / (2*pow(sigma, 2.0)));
}

double cdf_rayleigh(double x, double sigma) {
  if (sigma <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x == INFINITY)
    return 1;
  if (x >= 0)
    return 1 - exp(-pow(x, 2.0) / (2*pow(sigma, 2.0)));
  else
    return 0;
}

double invcdf_rayleigh(double p, double sigma) {
  if (p < 0 || p > 1 || sigma <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return sqrt(-2*pow(sigma, 2.0) * log(1-p));
}


// [[Rcpp::export]]
NumericVector cpp_drayleigh(
    const NumericVector& x,
    const NumericVector& sigma,
    bool log_prob = false
  ) {

  int n = x.length();
  int ns = sigma.length();
  int Nmax = std::max(n, ns);
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_rayleigh(x[i % n], sigma[i % ns]);

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_prayleigh(
    const NumericVector& x,
    const NumericVector& sigma,
    bool lower_tail = true, bool log_prob = false
  ) {

  int n = x.length();
  int ns = sigma.length();
  int Nmax = std::max(n, ns);
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_rayleigh(x[i % n], sigma[i % ns]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qrayleigh(
    const NumericVector& p,
    const NumericVector& sigma,
    bool lower_tail = true, bool log_prob = false
  ) {

  int n  = p.length();
  int ns = sigma.length();
  int Nmax = std::max(n, ns);
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);

  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1-pp[i];

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_rayleigh(pp[i % n], sigma[i % ns]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rrayleigh(
    const int n,
    const NumericVector& sigma
  ) {

  double u;
  int ns = sigma.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    u = R::runif(0, 1);
    x[i] = invcdf_rayleigh(u, sigma[i % ns]);
  }

  return x;
}

