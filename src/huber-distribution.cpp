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


double pdf_huber(double x, double mu, double sigma, double c) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(c))
    return NA_REAL;
  if (sigma <= 0.0 || c <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  double z, A, rho;
  z = abs((x - mu)/sigma);
  A = 2.0*SQRT_2_PI * (Phi(c) + phi(c)/c - 0.5);

  if (z <= c)
    rho = pow(z, 2.0)/2.0;
  else
    rho = c*z - pow(c, 2.0)/2.0;

  return exp(-rho)/A/sigma;
}

double cdf_huber(double x, double mu, double sigma, double c) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(c))
    return NA_REAL;
  if (sigma <= 0.0 || c <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }

  double A, z, az, p;
  A = 2.0*(phi(c)/c - Phi(-c) + 0.5);
  z = (x - mu)/sigma;
  az = -abs(z);
  
  if (az <= -c) 
    p = exp(pow(c, 2.0)/2.0)/c * exp(c*az) / SQRT_2_PI/A;
  else
    p = (phi(c)/c + Phi(az) - Phi(-c))/A;
  
  if (z <= 0.0)
    return p;
  else
    return 1.0 - p;
}

double invcdf_huber(double p, double mu, double sigma, double c) {
  if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma) || ISNAN(c))
    return NA_REAL;
  if (sigma <= 0.0 || c <= 0.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }

  double x, pm, A;
  A = 2.0 * SQRT_2_PI * (Phi(c) + phi(c)/c - 0.5);
  pm = std::min(p, 1.0 - p);

  if (pm <= SQRT_2_PI * phi(c)/(c*A))
    x = log(c*pm*A)/c - c/2.0;
  else
    x = InvPhi(abs(1.0 - Phi(c) + pm*A/SQRT_2_PI - phi(c)/c));

  if (p < 0.5)
    return mu + x*sigma;
  else
    return mu - x*sigma;
}

double rng_huber(double mu, double sigma, double c) {
  if (ISNAN(mu) || ISNAN(sigma) || ISNAN(c) || sigma <= 0.0 || c <= 0.0) {
    Rcpp::warning("NAs produced");
    return NA_REAL;
  }
  
  double x, pm, A, u;
  u = rng_unif();
  A = 2.0 * SQRT_2_PI * (Phi(c) + phi(c)/c - 0.5);
  pm = std::min(u, 1.0 - u);
  
  if (pm <= SQRT_2_PI * phi(c)/(c*A))
    x = log(c*pm*A)/c - c/2.0;
  else
    x = InvPhi(abs(1.0 - Phi(c) + pm*A/SQRT_2_PI - phi(c)/c));
  
  if (u < 0.5)
    return mu + x*sigma;
  else
    return mu - x*sigma;
}


// [[Rcpp::export]]
NumericVector cpp_dhuber(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& epsilon,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  dims.push_back(epsilon.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_huber(x[i % dims[0]], mu[i % dims[1]],
                     sigma[i % dims[2]], epsilon[i % dims[3]]);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_phuber(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& epsilon,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  dims.push_back(epsilon.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_huber(x[i % dims[0]], mu[i % dims[1]],
                     sigma[i % dims[2]], epsilon[i % dims[3]]);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qhuber(
    const NumericVector& p,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& epsilon,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  dims.push_back(epsilon.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_huber(pp[i % dims[0]], mu[i % dims[1]],
                        sigma[i % dims[2]], epsilon[i % dims[3]]);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rhuber(
    const int& n,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& epsilon
  ) {
  
  std::vector<int> dims;
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  dims.push_back(epsilon.length());
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_huber(mu[i % dims[0]], sigma[i % dims[1]], epsilon[i % dims[2]]);
  
  return x;
}

