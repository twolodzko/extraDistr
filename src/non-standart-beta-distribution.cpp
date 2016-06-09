#include <Rcpp.h>

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
*  Non-standard beta distribution
*
*  Values:
*  x
*
*  Parameters:
*  0 <= beta <= 1
*  alpha > 0
*  lower < upper
*
*/

double pdf_nsbeta(double x, double alpha, double beta, double l, double u, bool log_p) {
  if (l >= u || alpha < 0 || beta < 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double r = u-l;
  double p = R::dbeta((x-l)/r, alpha, beta, log_p);
  if (log_p) 
    return p-log(r);
  else
    return p/r;
}

double cdf_nsbeta(double x, double alpha, double beta, double l, double u, bool lower_tail, bool log_p) {
  if (l >= u || alpha < 0 || beta < 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return R::pbeta((x-l)/(u-l), alpha, beta, lower_tail, log_p);
}

double invcdf_nsbeta(double p, double alpha, double beta, double l, double u, bool lower_tail, bool log_p) {
  if (l >= u || alpha < 0 || beta < 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return R::qbeta(p, alpha, beta, lower_tail, log_p) * (u-l) + l;
}

double rng_nsbeta(double alpha, double beta, double l, double u) {
  if (l >= u || alpha < 0 || beta < 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return R::rbeta(alpha, beta) * (u-l) + l;
}


// [[Rcpp::export]]
NumericVector cpp_dnsbeta(
    const NumericVector& x,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& lower,
    const NumericVector& upper,
    bool log_prob = false
  ) {
  
  int n  = x.length();
  int na = beta.length();
  int nb = alpha.length();
  int nl = lower.length();
  int nu = upper.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb, nl, nu));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_nsbeta(x[i % n], alpha[i % na], beta[i % nb], lower[i % nl], upper[i % nu], log_prob);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pnsbeta(
    const NumericVector& x,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& lower,
    const NumericVector& upper,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = x.length();
  int na = beta.length();
  int nb = alpha.length();
  int nl = lower.length();
  int nu = upper.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb, nl, nu));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_nsbeta(x[i % n], alpha[i % na], beta[i % nb], lower[i % nl], upper[i % nu], lower_tail, log_prob);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qnsbeta(
    const NumericVector& p,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& lower,
    const NumericVector& upper,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = p.length();
  int na = beta.length();
  int nb = alpha.length();
  int nl = lower.length();
  int nu = upper.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb, nl, nu));
  NumericVector q(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_nsbeta(p[i % n], alpha[i % na], beta[i % nb], lower[i % nl], upper[i % nu], lower_tail, log_prob);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rnsbeta(
    const int n,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& lower,
    const NumericVector& upper
  ) {
  
  int na = beta.length();
  int nb = alpha.length();
  int nl = lower.length();
  int nu = upper.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_nsbeta(alpha[i % na], beta[i % nb], lower[i % nl], upper[i % nu]);
  
  return x;
}

