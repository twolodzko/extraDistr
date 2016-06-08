#include <Rcpp.h>
#include "namespace.h"
#include "const.h"
#include "shared.h"


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
  if (alpha <= 0 || beta <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(x) || x < 0 || std::isinf(x))
    return -INFINITY;
  
  double p = beta/(1+beta);
  return R::lgammafn(alpha+x) - (lfactorial(x) + R::lgammafn(alpha)) +
    log(p)*x + log(1-p)*alpha;
}

double cdf_gpois(double x, double alpha, double beta) {
  if (alpha <= 0 || beta <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (std::isinf(x))
    return 1;
  double p_tmp = 0;
  for (int j = 0; j < x+1; j++)
    p_tmp += exp(logpmf_gpois(j, alpha, beta))*P_NORM_CONST;
  return p_tmp/P_NORM_CONST;
}

double rng_gpois(double alpha, double beta) {
  if (alpha <= 0 || beta <= 0) {
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
    
    double mx = finite_max(x);
    NumericVector p_tab(mx+1);
    
    p_tab[0] = exp(logpmf_gpois(0, alpha[0], beta[0]))*P_NORM_CONST;
    for (int j = 1; j < mx+1; j++)
      p_tab[j] = p_tab[j-1] + exp(logpmf_gpois(j, alpha[0], beta[0]))*P_NORM_CONST;
    
    for (int i = 0; i < n; i++) {
      if (std::isinf(x[i])) {
        p[i] = 1;
      } else if (x[i] == floor(x[i]) && x[i] >= 0) {
        p[i] = p_tab[(int)x[i]]/P_NORM_CONST;
      } else {
        p[i] = 0;
      }
    }
    
  } else {
    
    for (int i = 0; i < Nmax; i++)
      p[i] = cdf_gpois(x[i % n], alpha[i % na], beta[i % nb]);
    
  }
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

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

