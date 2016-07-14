#include <Rcpp.h>
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
 * Birnbaum-Saunders (Fatigue Life) Distribution
 * 
 * Support:
 * x > mu
 * 
 * Parameters:
 * mu
 * alpha > 0
 * beta > 0
 * 
 * 
 */

double pdf_fatigue(double x, double alpha, double beta, double mu) {
  if (alpha <= 0.0 || beta <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x <= mu || std::isinf(x))
    return 0.0;
  double z, zb, bz;
  z = x-mu;
  zb = sqrt(z/beta);
  bz = sqrt(beta/z);
  return (zb+bz)/(2.0*alpha*z) * phi((zb-bz)/alpha);
}

double cdf_fatigue(double x, double alpha, double beta, double mu) {
  if (alpha <= 0.0 || beta <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x <= mu)
    return 0.0;
  double z, zb, bz;
  z = x-mu;
  zb = sqrt(z/beta);
  bz = sqrt(beta/z);
  return Phi((zb-bz)/alpha);
}

double invcdf_fatigue(double p, double alpha, double beta, double mu) {
  if (alpha <= 0.0 || beta <= 0.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (p == 0.0)
    return mu;
  double Zp = InvPhi(p);
  return pow(alpha/2.0*Zp + sqrt(pow(alpha/2.0*Zp, 2.0) + 1.0), 2.0) * beta + mu;
}

double rng_fatigue(double alpha, double beta, double mu) {
  if (alpha <= 0.0 || beta <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double z = R::norm_rand();
  return pow(alpha/2.0*z + sqrt(pow(alpha/2.0*z, 2.0) + 1.0), 2.0) * beta + mu;
}


// [[Rcpp::export]]
NumericVector cpp_dfatigue(
    const NumericVector& x,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& mu,
    bool log_prob = false
  ) {
  
  int n  = x.length();
  int na = alpha.length();
  int nb = beta.length();
  int nm = mu.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb, nm));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_fatigue(x[i % n], alpha[i % na], beta[i % nb], mu[i % nm]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pfatigue(
    const NumericVector& x,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& mu,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = x.length();
  int na = alpha.length();
  int nb = beta.length();
  int nm = mu.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb, nm));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_fatigue(x[i % n], alpha[i % na], beta[i % nb], mu[i % nm]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qfatigue(
    const NumericVector& p,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& mu,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = p.length();
  int na = alpha.length();
  int nb = beta.length();
  int nm = mu.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb, nm));
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1.0 - pp[i];
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_fatigue(pp[i % n], alpha[i % na], beta[i % nb], mu[i % nm]);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rfatigue(
    const int n,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& mu
  ) {
  
  int na = alpha.length();
  int nb = beta.length();
  int nm = mu.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_fatigue(alpha[i % na], beta[i % nb], mu[i % nm]);
  
  return x;
}

