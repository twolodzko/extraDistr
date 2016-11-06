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
  if (ISNAN(x) || ISNAN(alpha) || ISNAN(beta))
    return NA_REAL;
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

/*
double cdf_gpois(double x, double alpha, double beta) {
  if (ISNAN(x) || ISNAN(alpha) || ISNAN(beta))
    return NA_REAL;
  if (alpha <= 0.0 || beta <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x == INFINITY)
    return 1.0;
  if (x < 0.0)
    return 0.0;
  
  double p_tmp, p, qa, ga, gax, xf, px;
  
  p_tmp = 0.0;
  p = beta/(1.0+beta);
  qa = pow(1.0 - p, alpha);
  ga = R::gammafn(alpha);
  
  // x = 0
  
  gax = ga;
  xf = 1.0;
  px = 1.0;
  p_tmp += qa;
  
  if (x < 1.0)
    return p_tmp;
  
  // x < 2
  
  gax *= alpha;
  px *= p;
  p_tmp += gax/ga * px * qa;
  
  if (x < 2.0)
    return p_tmp;
  
  // x >= 2
  
  double i = 2.0;
  while (i <= x) {
    gax *= i + alpha - 1.0;
    xf *= i;
    px *= p;
    p_tmp += gax/(xf * ga) * px * qa;
    i += 1.0;
  }

  return p_tmp;
}
*/

std::vector<double> cdf_gpois_table(double x, double alpha, double beta) {
  
  x = floor(x);
  std::vector<double> p_tab(static_cast<int>(x)+1);
  double p, qa, ga, gax, xf, px, lp;
  
  p = beta/(1.0+beta);
  qa = log(pow(1.0 - p, alpha));
  ga = R::lgammafn(alpha);
  lp = log(p);
  
  // x = 0
  
  gax = ga;
  xf = 0.0;
  px = 0.0;
  p_tab[0] = exp(qa);
  
  if (x < 1.0)
    return p_tab;
  
  // x < 2
  
  gax += log(alpha);
  px += lp;
  p_tab[1] = p_tab[0] + exp(gax - ga + px + qa);
  
  if (x < 2.0)
    return p_tab;
  
  // x >= 2
  
  double i = 2.0;
  while (i <= x) {
    gax += log(i + alpha - 1.0);
    xf += log(i);
    px += lp;
    p_tab[static_cast<int>(i)] = p_tab[static_cast<int>(i)-1] + exp(gax - (xf + ga) + px + qa);
    i += 1.0;
  }
  
  return p_tab;
}

double rng_gpois(double alpha, double beta) {
  if (ISNAN(alpha) || ISNAN(beta))
    return NA_REAL;
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
  
  if (na == 1 && nb == 1) {
    
    if (ISNAN(alpha[0]) || ISNAN(beta[0]) || allNA(x)) {
      for (int i = 0; i < n; i++)
        p[i] = NA_REAL;
      return p;
    }
    
    if (alpha[0] <= 0.0 || beta[0] <= 0.0) {
      for (int i = 0; i < n; i++)
        p[i] = NAN;
      Rcpp::warning("NaNs produced");
      return p;
    }
    
    double mx = finite_max(x);
    if (mx < 0.0 || mx == INFINITY)
      mx = 0.0;
    std::vector<double> p_tab = cdf_gpois_table(mx, alpha[0], beta[0]);
    
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
    
    for (int i = 0; i < Nmax; i++) {
      if (i % 1000 == 0)
        Rcpp::checkUserInterrupt();
      if (ISNAN(x[i % n]) || ISNAN(alpha[i % na]) || ISNAN(beta[i % nb])) {
        p[i] = NA_REAL;
      } else if (alpha[i % na] <= 0.0 || beta[i % nb] <= 0.0) {
        Rcpp::warning("NaNs produced");
        p[i] = NAN;
      } else if (x[i % n] < 0.0) {
        p[i] = 0.0;
      } else if (x[i % n] == INFINITY) {
        p[i] = 1.0;
      } else {
        p[i] = cdf_gpois_table(x[i % n], alpha[i % na], beta[i % nb]).back();
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

