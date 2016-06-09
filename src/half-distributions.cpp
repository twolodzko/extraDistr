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
 * 
 * x >= 0
 * 
 * Parameters:
 * nu > 0
 * sigma > 0
 * 
 * with nu = 1   returns half-Cauchy
 * with nu = Inf returns half-normal
 * 
 */

// Half-Cauchy

double pdf_hcauchy(double x, double sigma) {
  if (sigma <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0)
    return 0;
  return 2/(M_PI*(1 + pow(x/sigma, 2.0)))/sigma;
}

double cdf_hcauchy(double x, double sigma) {
  if (sigma <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0)
    return 0;
  return 2/M_PI * atan(x/sigma);
}

double invcdf_hcauchy(double p, double sigma) {
  if (sigma <= 0 || p < 0 || p > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return sigma * tan((M_PI*p)/2);
}

double rng_hcauchy(double sigma) {
  if (sigma <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return abs(R::rcauchy(0, sigma));
}

// Half-t

double pdf_ht(double x, double nu, double sigma) {
  if (sigma <= 0 || nu <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0)
    return 0;
  return 2 * R::dt(x/sigma, nu, false)/sigma;
}

double cdf_ht(double x, double nu, double sigma) {
  if (sigma <= 0 || nu <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0)
    return 0;
  return 2 * R::pt(x/sigma, nu, true, false) - 1;
}

double invcdf_ht(double p, double nu, double sigma) {
  if (sigma <= 0 || nu <= 0 || p < 0 || p > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return R::qt((p+1)/2, nu, true, false) * sigma;
}

double rng_ht(double nu, double sigma) {
  if (sigma <= 0 || nu <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return abs(R::rt(nu) * sigma);
}

// Half-normal

double pdf_hnorm(double x, double sigma) {
  if (sigma <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0)
    return 0;
  return 2 * R::dnorm(x, 0, sigma, false);
}

double cdf_hnorm(double x, double sigma) {
  if (sigma <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x < 0)
    return 0;
  return 2 * R::pnorm(x, 0.0, sigma, true, false) - 1;
}

double invcdf_hnorm(double p, double sigma) {
  if (sigma <= 0 || p < 0 || p > 1) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return R::qnorm((p+1)/2, 0.0, sigma, true, false);
}

double rng_hnorm(double sigma) {
  if (sigma <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return abs(R::rnorm(0.0, sigma));
}


// [[Rcpp::export]]
NumericVector cpp_dhalf(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& sigma,
    bool log_prob = false
  ) {
  
  int n  = x.length();
  int nn = nu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, ns));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++) {
    if (nu[i % nn] == 1) {
      p[i] = pdf_hcauchy(x[i % n], sigma[i % ns]);
    } else if (nu[i % nn] == INFINITY) {
      p[i] = pdf_hnorm(x[i % n], sigma[i % ns]);
    } else {
      p[i] = pdf_ht(x[i % n], nu[i % nn], sigma[i % ns]);
    }
  }
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_phalf(
    const NumericVector& x,
    const NumericVector& nu,
    const NumericVector& sigma,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = x.length();
  int nn = nu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, ns));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++) {
    if (nu[i % nn] == 1) {
      p[i] = cdf_hcauchy(x[i % n], sigma[i % ns]);
    } else if (nu[i % nn] == INFINITY) {
      p[i] = cdf_hnorm(x[i % n], sigma[i % ns]);
    } else {
      p[i] = cdf_ht(x[i % n], nu[i % nn], sigma[i % ns]);
    }
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
NumericVector cpp_qhalf(
    const NumericVector& p,
    const NumericVector& nu,
    const NumericVector& sigma,
    bool lower_tail = true, bool log_prob = false
  ) {
  
  int n  = p.length();
  int nn = nu.length();
  int ns = sigma.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, ns));
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(pp[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1-pp[i];
  
  for (int i = 0; i < Nmax; i++) {
    if (nu[i % nn] == 1) {
      q[i] = invcdf_hcauchy(pp[i % n], sigma[i % ns]);
    } else if (nu[i % nn] == INFINITY) {
      q[i] = invcdf_hnorm(pp[i % n], sigma[i % ns]);
    } else {
      q[i] = invcdf_ht(pp[i % n], nu[i % nn], sigma[i % ns]);
    }
  }
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rhalf(
    const int n,
    const NumericVector& nu,
    const NumericVector& sigma
  ) {
  
  int nn = nu.length();
  int ns = sigma.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++) {
    if (nu[i % nn] == 1) {
      x[i] = rng_hcauchy(sigma[i % ns]);
    } else if (nu[i % nn] == INFINITY) {
      x[i] = rng_hnorm(sigma[i % ns]);
    } else {
      x[i] = rng_ht(nu[i % nn], sigma[i % ns]);
    }
  }
  
  return x;
}

