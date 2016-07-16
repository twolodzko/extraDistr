#include <Rcpp.h>
#include "const.h"
#include "shared.h"

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
*  Gamma-Poisson distribution
*
*  Values:
*  x >= 0
*
*  Parameters:
*  alpha > 0
*  beta > 0
*
*/

double logpmf_gpois(double x, double alpha, double beta) {
  if (alpha <= 0.0 || beta <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(x) || x < 0.0 || std::isinf(x))
    return -INFINITY;
  
  double p = beta/(1.0+beta);
  return R::lgammafn(alpha+x) - (lfactorial(x) + R::lgammafn(alpha)) +
    log(p)*x + log(1.0-p)*alpha;
}

double cdf_gpois(double x, double alpha, double beta) {
  if (alpha <= 0.0 || beta <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x == INFINITY)
    return 1.0;
  if (x < 0.0)
    return 0.0;
  double p_tmp = 0.0;
  for (int j = 0; j < static_cast<int>(x)+1; j++)
    p_tmp += exp(logpmf_gpois(static_cast<double>(j), alpha, beta));
  return p_tmp;
}

double rng_gpois(double alpha, double beta) {
  if (alpha <= 0.0 || beta <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double lambda = R::rgamma(alpha, beta);
  return R::rpois(lambda);
}


// [[Rcpp::export]]
NumericVector cpp_dgpois(
    const NumericVector& x,
    const NumericVector& alpha,
    const NumericVector& beta,
    bool log_prob = false
  ) {

  int n  = x.length();
  int na = alpha.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = logpmf_gpois(x[i % n], alpha[i % na], beta[i % nb]);

  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = exp(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pgpois(
    const NumericVector& x,
    const NumericVector& alpha,
    const NumericVector& beta,
    bool lower_tail = true, bool log_prob = false
  ) {

  int n  = x.length();
  int na = alpha.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector p(Nmax);

  if (na == 1 && nb == 1 && anyFinite(x)) {
    
    double mx = static_cast<int>(finite_max(x));
    NumericVector p_tab(mx+1);
    
    p_tab[0] = exp(logpmf_gpois(0, alpha[0], beta[0]));
    for (int j = 1; j < mx+1; j++)
      p_tab[j] = p_tab[j-1] + exp(logpmf_gpois(static_cast<double>(j),
                                               alpha[0], beta[0]));
    
    for (int i = 0; i < n; i++) {
      if (x[i] == INFINITY) {
        p[i] = 1.0;
      } else if (x[i] >= 0.0) {
        p[i] = p_tab[static_cast<int>(x[i])];
      } else {
        p[i] = 0.0;
      }
    }
    
  } else {
    
    for (int i = 0; i < Nmax; i++) {
      if (i % 1000 == 0)
        Rcpp::checkUserInterrupt();
      p[i] = cdf_gpois(x[i % n], alpha[i % na], beta[i % nb]);
    }
    
  }
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rgpois(
    const int n,
    const NumericVector& alpha,
    const NumericVector& beta
  ) {

  int na = alpha.length();
  int nb = beta.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_gpois(alpha[i % na], beta[i % nb]);

  return x;
}

