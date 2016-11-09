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
  if (ISNAN(x) || ISNAN(alpha) || ISNAN(beta) || ISNAN(mu))
    return NA_REAL;
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
  if (ISNAN(x) || ISNAN(alpha) || ISNAN(beta) || ISNAN(mu))
    return NA_REAL;
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
  if (ISNAN(p) || ISNAN(alpha) || ISNAN(beta) || ISNAN(mu))
    return NA_REAL;
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
  if (ISNAN(alpha) || ISNAN(beta) || ISNAN(mu))
    return NA_REAL;
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
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(alpha.length());
  dims.push_back(beta.length());
  dims.push_back(mu.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_fatigue(x[i % dims[0]], alpha[i % dims[1]],
                       beta[i % dims[2]], mu[i % dims[3]]);
  
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
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(alpha.length());
  dims.push_back(beta.length());
  dims.push_back(mu.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_fatigue(x[i % dims[0]], alpha[i % dims[1]],
                       beta[i % dims[2]], mu[i % dims[3]]);
  
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
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(alpha.length());
  dims.push_back(beta.length());
  dims.push_back(mu.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    for (int i = 0; i < dims[0]; i++)
      pp[i] = exp(pp[i]);
  
  if (!lower_tail)
    for (int i = 0; i < dims[0]; i++)
      pp[i] = 1.0 - pp[i];
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_fatigue(pp[i % dims[0]], alpha[i % dims[1]],
                          beta[i % dims[2]], mu[i % dims[3]]);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rfatigue(
    const int& n,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& mu
  ) {
  
  std::vector<int> dims;
  dims.push_back(alpha.length());
  dims.push_back(beta.length());
  dims.push_back(mu.length());
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_fatigue(alpha[i % dims[0]], beta[i % dims[1]], mu[i % dims[2]]);
  
  return x;
}

