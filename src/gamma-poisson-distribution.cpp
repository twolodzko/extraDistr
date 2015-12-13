#include <Rcpp.h>
#include "shared.h"
using namespace Rcpp;


/*
*  Gamma-Poisson distribution
*
*  Values:
*  x >= 0
*
*  Parameters:
*  alpha > 0
*  beta > 0
*
*  f(x) = (gamma(x+beta)*alpha^x) / (gamma(beta)*(1+alpha)^(beta+x) * x!)
*  
*
*/

// double pmf_gpois(double x, double alpha, double beta) {
//   if (std::floor(x) != x)
//     return 0;
//   if (alpha <= 0 || beta <= 0)
//     return NAN;
//   if (x >= 0) {
//     return (R::gammafn(x+beta) * std::pow(alpha, x)) /
//            (R::gammafn(beta) * std::pow(1+alpha, beta+x) * factorial(x));
//   } else {
//     return 0;
//   }
// }

double logpmf_gpois2(double x, double alpha, double beta) {
  if (std::floor(x) != x)
    return 0;
  if (alpha <= 0 || beta <= 0)
    return NAN;
  if (x >= 0) {
    double p = beta/(1+beta);
    return R::lgammafn(alpha+x) - (lfactorial(x) + R::lgammafn(alpha)) +
           std::log(p)*x + std::log(1-p)*alpha;
  } else {
    return 0;
  }
}

// double cdf_gpois(double x, double alpha, double beta) {
//   if (alpha <= 0 || beta <= 0)
//     return NAN;
//   double p_tmp = 0;
//   for (int j = 0; j < x+1; j++)
//     p_tmp += pmf_gpois(j, alpha, beta);
//   return p_tmp;
// }

double rng_gpois(double alpha, double beta) {
  double lambda = R::rgamma(alpha, beta);
  return R::rpois(lambda);
}

// double logpmf_gpois(double x, double alpha, double beta) {
//   if (std::floor(x) != x)
//     return -INFINITY;
//   if (alpha <= 0 || beta <= 0)
//     return NAN;
//   if (x >= 0) {
//     return (R::lgammafn(x+beta) + std::log(alpha)*x) -
//            (R::lgammafn(beta) + std::log(1+alpha)*(beta+x) + lfactorial(x));
//   } else {
//     return -INFINITY;
//   }
// }
// 
// double cdf2_gpois(int x, double alpha, double beta) {
//   if (alpha <= 0 || beta <= 0)
//     return NAN;
//   double p_tmp = 0;
//   for (int j = 0; j < x+1; j++)
//     p_tmp += std::exp(logpmf_gpois(j, alpha, beta));
//   return p_tmp;
// }


// [[Rcpp::export]]
NumericVector cpp_dgpois(NumericVector x,
                         NumericVector alpha, NumericVector beta,
                         bool log_prob = false) {

  int n  = x.length();
  int na = alpha.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = logpmf_gpois2(x[i % n], alpha[i % na], beta[i % nb]);

  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::exp(p[i]);

  return p;
}

// 
// // [[Rcpp::export]]
// NumericVector cpp_pgpois(NumericVector x,
//                          NumericVector alpha, NumericVector beta,
//                          bool lower_tail = true, bool log_prob = false) {
// 
//   int n  = x.length();
//   int na = alpha.length();
//   int nb = beta.length();
//   int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
//   NumericVector p(Nmax);
// 
//   for (int i = 0; i < Nmax; i++)
//     p[i] = cdf2_gpois(x[i % n], alpha[i % na], beta[i % nb]);
// 
//   if (!lower_tail)
//     for (int i = 0; i < Nmax; i++)
//       p[i] = 1-p[i];
// 
//   if (log_prob)
//     for (int i = 0; i < Nmax; i++)
//       p[i] = std::log(p[i]);
// 
//   return p;
// }


// [[Rcpp::export]]
NumericVector cpp_rgpois(int n,
                         NumericVector alpha, NumericVector beta) {

  int na = alpha.length();
  int nb = beta.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_gpois(alpha[i % na], beta[i % nb]);

  return x;
}

