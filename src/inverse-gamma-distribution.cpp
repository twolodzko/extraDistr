#include <Rcpp.h>
using namespace Rcpp;


/*
*  Inverse-Gamma distribution
*
*  Values:
*  x
*
*  Parameters:
*  alpha > 0
*  beta > 0
*
*  f(k) = (x^(-alpha-1) * exp(-1/(beta*x))) / (Gamma(alpha)*beta^alpha)
*  F(x) = Gamma(alpha, 1/(beta*x)) / Gamma(alpha)
*
*  V. Witkovsky (2001) Computing the distribution of a linear
*  combination of inverted gamma variables, Kybernetika 37(1), 79-90
*
*/

double pdf_invgamma(double x, double alpha, double beta) {
  if (x > 0)
    return (pow(x, -alpha-1) * exp(-1/(beta*x))) / (R::gammafn(alpha) * pow(beta, alpha));
  else
    return 0;
}

double cdf_invgamma(double x, double alpha, double beta) {
  if (x > 0)
    return (R::pgamma(1/(beta*x), alpha, 1, false, false) * R::gammafn(alpha)) / R::gammafn(alpha);
  else
    return 0;
}

double rng_invgamma(double alpha, double beta) {
  return 1/R::rgamma(alpha, beta);
}

double logcdf_invgamma(double x, double alpha, double beta) {
  if (x > 0)
    return (R::pgamma(1/(beta*x), alpha, 1, false, true) + R::lgammafn(alpha)) - R::lgammafn(alpha);
  else
    return -INFINITY;
}


// [[Rcpp::export]]
NumericVector cpp_dinvgamma(NumericVector x,
                            NumericVector alpha, NumericVector beta,
                            bool log_prob = false) {

  if (is_true(any(alpha <= 0)) || is_true(any(beta <= 0)))
    throw Rcpp::exception("Values of alpha and beta should be > 0.");

  int n = x.length();
  int na = alpha.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_invgamma(x[i % n], alpha[i % na], beta[i % nb]);

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pinvgamma(NumericVector x,
                            NumericVector alpha, NumericVector beta,
                            bool lower_tail = true, bool log_prob = false) {

  if (is_true(any(alpha <= 0)) || is_true(any(beta <= 0)))
    throw Rcpp::exception("Values of alpha and beta should be > 0.");

  int n = x.length();
  int na = alpha.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_invgamma(x[i % n], alpha[i % na], beta[i % nb]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rinvgamma(int n,
                            NumericVector alpha, NumericVector beta) {

  if (is_true(any(alpha <= 0)) || is_true(any(beta <= 0)))
    throw Rcpp::exception("Values of alpha and beta should be > 0.");

  int na = alpha.length();
  int nb = beta.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_invgamma(alpha[i % na], beta[i % nb]);

  return x;
}

