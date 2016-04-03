#include <Rcpp.h>
#include "const.h"
#include "shared.h"
using namespace Rcpp;

/*
*  Beta-negative binomial distribution
*
*  Values:
*  x
*
*  Parameters:
*  r > 0
*  alpha > 0
*  beta > 0
*
*  f(k) = gamma(r+k)/(k! gamma(r)) * beta(alpha+r, beta+k)/beta(alpha, beta)
*
*/

double pmf_bnbinom(double k, double r, double alpha, double beta) {
  if (alpha <= 0 || beta <= 0 || floor(r) != r) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(k))
    return 0;
  return (R::gammafn(r+k) / (R::gammafn(k+1) * R::gammafn(r))) *
          R::beta(alpha+r, beta+k) / R::beta(alpha, beta);
}

double logpmf_bnbinom(double k, double r, double alpha, double beta) {
  if (alpha <= 0 || beta <= 0 || floor(r) != r) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(k))
    return -INFINITY;
  return (R::lgammafn(r+k) - (R::lgammafn(k+1) + R::lgammafn(r))) +
    R::lbeta(alpha+r, beta+k) - R::lbeta(alpha, beta);
}

double cdf_bnbinom(double k, double r, double alpha, double beta) {
  if (alpha < 0 || beta < 0 || floor(r) != r) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(k))
    return 0;
  double p_tmp = 0;
  for (int j = 0; j < k+1; j++)
    p_tmp += exp(logpmf_bnbinom(j, r, alpha, beta))*P_NORM_CONST;
  return p_tmp/P_NORM_CONST;
}

double rng_bnbinom(double r, double alpha, double beta) {
  if (alpha <= 0 || beta <= 0 || floor(r) != r) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double prob = R::rbeta(alpha, beta);
  return R::rnbinom(r, prob);
}


// [[Rcpp::export]]
NumericVector cpp_dbnbinom(NumericVector x,
                           NumericVector size, NumericVector alpha, NumericVector beta,
                           bool log_prob = false) {

  int n = x.length();
  int nn = size.length();
  int na = alpha.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = logpmf_bnbinom(x[i % n], size[i % nn], alpha[i % na], beta[i % nb]);

  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = exp(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pbnbinom(NumericVector x,
                           NumericVector size, NumericVector alpha, NumericVector beta,
                           bool lower_tail = true, bool log_prob = false) {

  int n = x.length();
  int nn = size.length();
  int na = alpha.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, na, nb));
  NumericVector p(Nmax);

  if (nn == 1 && na == 1 && nb == 1) {
    
    if (alpha[0] < 0 || beta[0] < 0 || floor(size[0]) != size[0]) {
      Rcpp::warning("NaNs produced");
      for (int i = 0; i < n; i++)
        p[i] = NAN;
      return p;
    }
    
    double mx = (int)floor(max(x));
    NumericVector p_tab(mx+1);
    
    p_tab[0] = exp(logpmf_bnbinom(0, size[0], alpha[0], beta[0]))*P_NORM_CONST;
    for (int j = 1; j < mx+1; j++)
      p_tab[j] = p_tab[j-1] + exp(logpmf_bnbinom(j, size[0], alpha[0], beta[0]))*P_NORM_CONST;
    
    for (int i = 0; i < n; i++) {
      if (x[i] == floor(x[i]) && x[i] >= 0)
        p[i] = p_tab[(int)x[i]]/P_NORM_CONST;
      else
        p[i] = 0;
    }
    
  } else {
    
    for (int i = 0; i < Nmax; i++)
      p[i] = cdf_bnbinom(x[i % n], size[i % nn], alpha[i % na], beta[i % nb]);
    
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
NumericVector cpp_rbnbinom(int n,
                           NumericVector size, NumericVector alpha, NumericVector beta) {

  int nn = size.length();
  int na = alpha.length();
  int nb = beta.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_bnbinom(size[i % nn], alpha[i % na], beta[i % nb]);

  return x;
}

