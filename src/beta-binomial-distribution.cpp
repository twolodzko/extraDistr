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
*  Beta-binomial distribution
*
*  Values:
*  x
*
*  Parameters:
*  k > 0
*  alpha > 0
*  beta > 0
*
*  f(k) = choose(n, k) * (beta(k+alpha, n-k+beta)) / (beta(alpha, beta))
*
*/

double pmf_bbinom(double k, double n, double alpha, double beta) {
  if (ISNAN(k) || ISNAN(n) || ISNAN(alpha) || ISNAN(beta))
    return NA_REAL;
  if (alpha < 0.0 || beta < 0.0 || n < 0.0 || floor(n) != n) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(k) || k < 0.0 || k > n)
    return 0.0;
  return R::choose(n, k) * R::beta(k+alpha, n-k+beta) / R::beta(alpha, beta);
}

double logpmf_bbinom(double k, double n, double alpha, double beta) {
  if (ISNAN(k) || ISNAN(n) || ISNAN(alpha) || ISNAN(beta))
    return NA_REAL;
  if (alpha < 0.0 || beta < 0.0 || n < 0.0 || floor(n) != n) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(k) || k < 0.0 || k > n)
    return -INFINITY;
  return R::lchoose(n, k) + R::lbeta(k+alpha, n-k+beta) - R::lbeta(alpha, beta);
}

/*
double cdf_bbinom(double k, double n, double alpha, double beta) {
  if (ISNAN(k) || ISNAN(n) || ISNAN(alpha) || ISNAN(beta))
    return NA_REAL;
  if (alpha < 0.0 || beta < 0.0 || n < 0.0 || floor(n) != n) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (k < 0.0)
    return 0.0;
  if (k > n)
    return 1.0;
  
  double p_tmp, nck, bab, gx, gy, gxy;
  
  p_tmp = 0.0;
  bab = R::beta(alpha, beta);
  gxy = R::gammafn(alpha + beta + n);

  // k = 0
  
  nck = 1.0;
  gx = R::gammafn(alpha);
  gy = R::gammafn(beta + n);
  p_tmp += nck * (gx * gy)/gxy/bab;
  
  if (k < 1.0)
    return p_tmp;
  
  // k < 2
  
  nck *= n;
  gx *= alpha;
  gy /= n + beta - 1.0;
  p_tmp += nck * (gx * gy)/gxy/bab;
  
  if (k < 2.0)
    return p_tmp;
  
  // k >= 1
  
  double i = 2.0;
  while (i <= k) {
    nck *= (n + 1.0 - i)/i;
    gx *= i + alpha - 1.0;
    gy /= n + beta - i;
    p_tmp += nck * (gx * gy)/gxy/bab;
    i += 1.0;
  }
  
  return p_tmp;
}
*/

std::vector<double> cdf_bbinom_table(double k, double n, double alpha, double beta) {

  k = floor(k);
  std::vector<double> p_tab(static_cast<int>(k)+1);
  double nck, bab, gx, gy, gxy;
  
  bab = R::lbeta(alpha, beta);
  gxy = R::lgammafn(alpha + beta + n);
  
  // k = 0
  
  nck = 0.0;
  gx = R::lgammafn(alpha);
  gy = R::lgammafn(beta + n);
  p_tab[0] = exp(nck + gx + gy - gxy - bab);
  
  if (k < 1.0)
    return p_tab;
  
  // k < 2
  
  nck += log(n);
  gx += log(alpha);
  gy -= log(n + beta - 1.0);
  p_tab[1] = p_tab[0] + exp(nck + gx + gy - gxy + bab);
  
  if (k < 2.0)
    return p_tab;
  
  // k >= 1
  
  double i = 2.0;
  while (i <= k) {
    nck += log((n + 1.0 - i)/i);
    gx += log(i + alpha - 1.0);
    gy -= log(n + beta - i);
    p_tab[static_cast<int>(i)] = p_tab[static_cast<int>(i)-1] + exp(nck + gx + gy - gxy - bab);
    i += 1.0;
  }
  
  return p_tab;
}

double rng_bbinom(double n, double alpha, double beta) {
  if (ISNAN(n) || ISNAN(alpha) || ISNAN(beta))
    return NA_REAL;
  if (alpha < 0.0 || beta < 0.0 || n < 0.0 || floor(n) != n) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double prob = R::rbeta(alpha, beta);
  return R::rbinom(n, prob);
}


// [[Rcpp::export]]
NumericVector cpp_dbbinom(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta,
    bool log_prob = false
  ) {

  int n = x.length();
  int nn = size.length();
  int na = alpha.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, na, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = logpmf_bbinom(x[i % n], size[i % nn], alpha[i % na], beta[i % nb]);

  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = exp(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pbbinom(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta,
    bool lower_tail = true, bool log_prob = false
  ) {

  int n = x.length();
  int nn = size.length();
  int na = alpha.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, na, nb));
  NumericVector p(Nmax);
  
  if (na == 1 && nb == 1 && nn == 1) {
    
    if (ISNAN(alpha[0]) || ISNAN(beta[0]) || allNA(x)) {
      for (int i = 0; i < n; i++)
        p[i] = NA_REAL;
      return p;
    }
    
    if (alpha[0] <= 0.0 || beta[0] <= 0.0  ||
        size[0] < 0.0 || floor(size[0]) != size[0]) {
      for (int i = 0; i < n; i++)
        p[i] = NAN;
      Rcpp::warning("NaNs produced");
      return p;
    }
    
    double mx = floor(finite_max(x));
    mx = std::max(0.0, std::min(mx, floor(size[0])));
    std::vector<double> p_tab = cdf_bbinom_table(mx, size[0], alpha[0], beta[0]);
    
    for (int i = 0; i < n; i++) {
      if (x[i] < 0.0) {
        p[i] = 0.0;
      } else if (x[i] >= size[0]) {
        p[i] = 1.0;
      } else {
        p[i] = p_tab[static_cast<int>(floor(x[i]))];
      } 
    }
    
  } else {
  
    double xi;
    for (int i = 0; i < Nmax; i++) {
      if (i % 1000 == 0)
        Rcpp::checkUserInterrupt();
      xi = floor(x[i % n]);
      if (ISNAN(xi) || ISNAN(size[i % nn]) || ISNAN(alpha[i % na]) || ISNAN(beta[i % nb])) {
        p[i] = NA_REAL;
      } else if (alpha[i % na] <= 0.0 || beta[i % nb] <= 0.0 ||
        size[i % nn] < 0.0 || floor(size[i % nn]) != size[i % nn]) {
        Rcpp::warning("NaNs produced");
        p[i] = NAN;
      } else if (xi < 0.0) {
        p[i] = 0.0;
      } else if (xi >= size[i % nn]) {
        p[i] = 1.0;
      } else {
        p[i] = cdf_bbinom_table(xi, size[i % nn], alpha[i % na], beta[i % nb]).back();
      }
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
NumericVector cpp_rbbinom(
    const int n,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta
  ) {

  int nn = size.length();
  int na = alpha.length();
  int nb = beta.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_bbinom(size[i % nn], alpha[i % na], beta[i % nb]);

  return x;
}

