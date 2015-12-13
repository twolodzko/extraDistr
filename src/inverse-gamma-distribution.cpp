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
*  F(x) = gamma(alpha, 1/(beta*x)) / Gamma(alpha)
*
*  V. Witkovsky (2001) Computing the distribution of a linear
*  combination of inverted gamma variables, Kybernetika 37(1), 79-90
*
*/

double pdf_invgamma(double x, double alpha, double beta) {
  if (alpha <= 0 || beta <= 0)
    return NAN;
  if (x > 0)
    return (std::pow(x, -alpha-1) * std::exp(-1/(beta*x))) / (R::gammafn(alpha) * std::pow(beta, alpha));
  else
    return 0;
}


// [[Rcpp::export]]
NumericVector cpp_dinvgamma(NumericVector x,
                            NumericVector alpha, NumericVector beta,
                            bool log_prob = false) {

  int n = x.length();
  int na = alpha.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_invgamma(x[i % n], alpha[i % na], beta[i % nb]);

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);

  return p;
}

