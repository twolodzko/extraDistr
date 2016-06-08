#include <Rcpp.h>
#include "namespace.h"
#include "const.h"
#include "shared.h"


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
    alpha_tot = 0;
    p[i] = 0;
    for (int j = 0; j < k; j++) {
      p[i] += alpha(i % na, j) * R::dnorm(x[i], mu(i % nm, j), sigma(i % ns, j), false);
      alpha_tot += alpha(i % na, j)*P_NORM_CONST;
      if (sigma(i % ns, j) <= 0)
        wrong_param = true;
    }
    if (!tol_equal(alpha_tot/P_NORM_CONST, 1) || wrong_param) {
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
    alpha_tot = 0;
    p[i] = 0;
    for (int j = 0; j < k; j++) {
      p[i] += alpha(i % na, j) * R::pnorm(x[i], mu(i % nm, j), sigma(i % ns, j), lower_tail, false);
      alpha_tot += alpha(i % na, j)*P_NORM_CONST;
      if (sigma(i % ns, j) <= 0)
        wrong_param = true;
    }
    if (!tol_equal(alpha_tot/P_NORM_CONST, 1) || wrong_param) {
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
  
  bool wrong_param, j_set;
  double u, p_tmp;
  int jj;
  NumericVector prob(k);
  
  for (int i = 0; i < n; i++) {
    wrong_param = false;
    j_set = false;
    u = R::runif(0, 1)*P_NORM_CONST;
    p_tmp = P_NORM_CONST;
    for (int j = k-1; j >= 0; j--) {
      p_tmp -= alpha(i % na, j)*P_NORM_CONST;
      if (sigma(i % ns, j) <= 0)
        wrong_param = true;
      if (!j_set && u > p_tmp) {
        jj = j;
        j_set = true;
      }
    }
    
    if (!tol_equal(p_tmp/P_NORM_CONST, 0) || wrong_param) {
      Rcpp::warning("NaNs produced");
      x[i] = NAN;
    } else {
      x[i] = R::rnorm(mu(i % nm, jj), sigma(i % ns, jj)); 
    }
  }
  
  return x;
}

