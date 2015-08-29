#include <Rcpp.h>
using namespace Rcpp;

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
  if (x >= 0)
    return x/pow(sigma, 2) * exp(-pow(x, 2) / (2*pow(sigma, 2)));
  else
    return 0;
}

double cdf_rayleigh(double x, double sigma) {
  if (x >= 0)
    return 1 - exp(-pow(x, 2) / pow(2*sigma, 2));
  else
    return 0;
}

double invcdf_rayleigh(double p, double sigma) {
  return sigma * sqrt(-2*log(1-p));
}


// [[Rcpp::export]]
NumericVector cpp_drayleigh(NumericVector x,
                            NumericVector sigma,
                            bool log_prob = false) {

  if (is_true(any(sigma <= 0)))
    throw Rcpp::exception("Values of sigma should be > 0.");

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
NumericVector cpp_prayleigh(NumericVector x,
                            NumericVector sigma,
                            bool lower_tail = true, bool log_prob = false) {

  if (is_true(any(sigma <= 0)))
    throw Rcpp::exception("Values of sigma should be > 0.");

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
NumericVector cpp_qrayleigh(NumericVector p,
                            NumericVector sigma,
                            bool lower_tail = true, bool log_prob = false) {

  if (is_true(any(p < 0)) || is_true(any(p > 1)))
    throw Rcpp::exception("Probabilities should range from 0 to 1.");
  if (is_true(any(sigma <= 0)))
    throw Rcpp::exception("Values of sigma should be > 0.");

  int n  = p.length();
  int ns = sigma.length();
  int Nmax = std::max(n, ns);
  NumericVector q(Nmax);

  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = exp(p[i]);

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      p[i] = 1-p[i];

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_rayleigh(p[i % n], sigma[i % ns]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rrayleigh(int n,
                            NumericVector sigma) {

  if (is_true(any(sigma <= 0)))
    throw Rcpp::exception("Values of sigma should be > 0.");

  double u;
  int ns = sigma.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    u = R::runif(0, 1);
    x[i] = invcdf_rayleigh(u, sigma[i % ns]);
  }

  return x;
}

