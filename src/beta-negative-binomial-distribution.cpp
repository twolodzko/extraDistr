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
  if (ISNAN(k) || ISNAN(r) || ISNAN(alpha) || ISNAN(beta))
    return NA_REAL;
  if (alpha <= 0.0 || beta <= 0.0 || r < 0.0 || floor(r) != r) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(k) || k < 0.0 || std::isinf(k))
    return 0.0;
  return (R::gammafn(r+k) / (R::gammafn(k+1.0) * R::gammafn(r))) *
          R::beta(alpha+r, beta+k) / R::beta(alpha, beta);
}

double logpmf_bnbinom(double k, double r, double alpha, double beta) {
  if (ISNAN(k) || ISNAN(r) || ISNAN(alpha) || ISNAN(beta))
    return NA_REAL;
  if (alpha <= 0.0 || beta <= 0.0 || r < 0.0 || floor(r) != r) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(k) || k < 0.0 || std::isinf(k))
    return -INFINITY;
  return (R::lgammafn(r+k) - (R::lgammafn(k+1.0) + R::lgammafn(r))) +
    R::lbeta(alpha+r, beta+k) - R::lbeta(alpha, beta);
}

/*
double cdf_bnbinom(double k, double r, double alpha, double beta) {
  if (ISNAN(k) || ISNAN(r) || ISNAN(alpha) || ISNAN(beta))
    return NA_REAL;
  if (alpha < 0.0 || beta < 0.0 || r < 0.0 || floor(r) != r) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (k < 0.0)
    return 0.0;
  if (k == INFINITY)
    return 1.0;

  double p_tmp, grx, xf, gr, gar, gbx, gabrx, bab;
  
  p_tmp = 0.0;
  bab = R::beta(alpha, beta);
  gr = R::gammafn(r);
  gar = R::gammafn(alpha + r);
  
  // k < 1
  
  grx = gr;
  gbx = R::gammafn(beta);
  gabrx = R::gammafn(alpha + beta + r);
  p_tmp += grx/gr * (gar*gbx)/gabrx/bab;
  
  if (k < 1.0)
    return p_tmp;
  
  // k < 2
  
  grx *= r;
  gbx *= beta;
  gabrx *= alpha + beta + r;
  p_tmp += grx/gr * (gar*gbx)/gabrx/bab;
  
  if (k < 2.0)
    return p_tmp;
  
  // k >= 2
  
  double i = 2.0;
  while (i <= k) {
    grx *= r + i - 1;
    gbx *= beta + i - 1;
    gabrx *= alpha + beta + r + i - 1;
    xf *= i;
    p_tmp += grx/(xf * gr) * (gar*gbx)/gabrx/bab;
    i += 1.0;
  }
  
  return p_tmp;
}
*/

std::vector<double> cdf_bnbinom_table(double k, double r, double alpha, double beta) {

  k = floor(k);
  std::vector<double> p_tab(static_cast<int>(k)+1);
  double grx, xf, gr, gar, gbx, gabrx, bab;
  
  bab = R::lbeta(alpha, beta);
  gr = R::lgammafn(r);
  gar = R::lgammafn(alpha + r);
  xf = 0.0;
  
  // k < 1
  
  grx = gr;
  gbx = R::lgammafn(beta);
  gabrx = R::lgammafn(alpha + beta + r);
  p_tab[0] = exp(grx - gr + gar + gbx - gabrx - bab);
  
  if (k < 1.0)
    return p_tab;
  
  // k < 2
  
  grx += log(r);
  gbx += log(beta);
  gabrx += log(alpha + beta + r);
  p_tab[1] = p_tab[0] + exp(grx - gr + gar + gbx - gabrx - bab);
  
  if (k < 2.0)
    return p_tab;
  
  // k >= 2
  
  double i = 2.0;
  while (i <= k) {
    grx += log(r + i - 1.0);
    gbx += log(beta + i - 1.0);
    gabrx += log(alpha + beta + r + i - 1.0);
    xf += log(i);
    p_tab[static_cast<int>(i)] = p_tab[static_cast<int>(i)-1] + exp(grx - (xf + gr) + gar + gbx - gabrx - bab);
    i += 1.0;
  }
  
  return p_tab;
}

double rng_bnbinom(double r, double alpha, double beta) {
  if (ISNAN(r) || ISNAN(alpha) || ISNAN(beta))
    return NA_REAL;
  if (alpha <= 0.0 || beta <= 0.0 || r < 0.0 || floor(r) != r) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double prob = R::rbeta(alpha, beta);
  return R::rnbinom(r, prob);
}


// [[Rcpp::export]]
NumericVector cpp_dbnbinom(
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
    p[i] = logpmf_bnbinom(x[i % n], size[i % nn], alpha[i % na], beta[i % nb]);

  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = exp(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pbnbinom(
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
    
    if (ISNAN(alpha[0]) || ISNAN(beta[0]) || ISNAN(size[0]) || allNA(x)) {
      for (int i = 0; i < n; i++)
        p[i] = NA_REAL;
      return p;
    }
    
    if (alpha[0] <= 0.0 || beta[0] <= 0.0 ||
        size[0] < 0.0 || floor(size[0]) != size[0]) {
      for (int i = 0; i < n; i++)
        p[i] = NAN;
      Rcpp::warning("NaNs produced");
      return p;
    }
    
    double mx = floor(finite_max(x));
    if (mx < 0.0 || mx == INFINITY)
      mx = 0.0;
    std::vector<double> p_tab = cdf_bnbinom_table(mx, size[0], alpha[0], beta[0]);
    
    for (int i = 0; i < n; i++) {
      if (ISNAN(x[i])) {
        p[i] = NA_REAL;
      } else if (x[i] < 0.0) {
        p[i] = 0.0;
      } else if (x[i] == INFINITY) {
        p[i] = 1.0;
      } else {
        p[i] = p_tab[static_cast<int>(x[i])];
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
      } else if (xi == INFINITY) {
        p[i] = 1.0;
      } else {
        p[i] = cdf_bnbinom_table(xi, size[i % nn], alpha[i % na], beta[i % nb]).back();
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
NumericVector cpp_rbnbinom(
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
    x[i] = rng_bnbinom(size[i % nn], alpha[i % na], beta[i % nb]);

  return x;
}

