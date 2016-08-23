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


// [[Rcpp::export]]
NumericVector cpp_dmixnorm(
    const NumericVector& x,
    const NumericMatrix& mu,
    const NumericMatrix& sigma,
    const NumericMatrix& alpha,
    bool log_prob = false
) {
  
  int n  = x.length();
  int nm = mu.nrow();
  int ns = sigma.nrow();
  int na = alpha.nrow();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, na));
  int k = alpha.ncol();
  NumericVector p(Nmax);
  
  if (k != mu.ncol() || k != sigma.ncol())
    Rcpp::stop("sizes of 'mu', 'sigma', and 'alpha' do not match");
  
  bool wrong_param;
  double alpha_tot;
  
  for (int i = 0; i < Nmax; i++) {
    wrong_param = false;
    alpha_tot = 0.0;
    p[i] = 0.0;
    
    for (int j = 0; j < k; j++) {
      if (alpha(i % na, j) < 0.0 || sigma(i % ns, j) <= 0.0) {
        wrong_param = true;
        break;
      }
      alpha_tot += alpha(i % na, j);
    }
    
    for (int j = 0; j < k; j++)
      p[i] += (alpha(i % na, j) / alpha_tot) * R::dnorm(x[i], mu(i % nm, j), sigma(i % ns, j), false);
    
    if (wrong_param) {
      Rcpp::warning("NaNs produced");
      p[i] = NAN;
    }
  }
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pmixnorm(
    const NumericVector& x,
    const NumericMatrix& mu,
    const NumericMatrix& sigma,
    const NumericMatrix& alpha,
    bool lower_tail = true, bool log_prob = false
) {
  
  int n  = x.length();
  int nm = mu.nrow();
  int ns = sigma.nrow();
  int na = alpha.nrow();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, na));
  int k = alpha.ncol();
  NumericVector p(Nmax);
  
  if (k != mu.ncol() || k != sigma.ncol())
    Rcpp::stop("sizes of 'mu', 'sigma', and 'alpha' do not match");
  
  bool wrong_param;
  double alpha_tot;
  
  for (int i = 0; i < Nmax; i++) {
    wrong_param = false;
    alpha_tot = 0.0;
    p[i] = 0.0;
    
    for (int j = 0; j < k; j++) {
      if (alpha(i % na, j) < 0.0 || sigma(i % ns, j) < 0.0) {
        wrong_param = true;
        break;
      }
      alpha_tot += alpha(i % na, j);
    }
    
    for (int j = 0; j < k; j++)
      p[i] += (alpha(i % na, j) / alpha_tot) * R::pnorm(x[i], mu(i % nm, j), sigma(i % ns, j), lower_tail, false);
    
    if (wrong_param) {
      Rcpp::warning("NaNs produced");
      p[i] = NAN;
    }
  }
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rmixnorm(
    const int n,
    const NumericMatrix& mu,
    const NumericMatrix& sigma,
    const NumericMatrix& alpha
) {
  
  int nm = mu.nrow();
  int ns = sigma.nrow();
  int na = alpha.nrow();
  int k = alpha.ncol();
  NumericVector x(n);
  
  if (k != mu.ncol() || k != sigma.ncol())
    Rcpp::stop("sizes of 'mu', 'sigma', and 'alpha' do not match");
  
  int jj;
  bool wrong_param;
  double u, p_tmp, alpha_tot;
  NumericVector prob(k);
  
  for (int i = 0; i < n; i++) {
    jj = 0;
    wrong_param = false;
    u = rng_unif();
    p_tmp = 1.0;
    alpha_tot = 0.0;
    
    for (int j = 0; j < k; j++) {
      if (alpha(i % na, j) < 0.0 || sigma(i % ns, j) < 0.0) {
        wrong_param = true;
        break;
      }
      alpha_tot += alpha(i % na, j);
    }
    
    for (int j = k-1; j >= 0; j--) {
      p_tmp -= alpha(i % na, j) / alpha_tot;
      if (u > p_tmp) {
        jj = j;
        break;
      }
    }

    if (wrong_param) {
      Rcpp::warning("NaNs produced");
      x[i] = NAN;
    } else {
      x[i] = R::rnorm(mu(i % nm, jj), sigma(i % ns, jj)); 
    }
  }
  
  return x;
}

