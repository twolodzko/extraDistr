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
*  Truncated Normal distribution
*
*  Values:
*  x
*
*  Parameters:
*  mu
*  sigma > 0
*  a, b
*
*  z = (x-mu)/sigma
*
*  f(x)    = phi(z) / (Phi((b-mu)/sigma) - Phi((mu-a)/sigma))
*  F(x)    = (Phi(z) - Phi((mu-a)/sigma)) / (Phi((b-mu)/sigma) - Phi((a-mu)/sigma))
*  F^-1(p) = Phi^-1(Phi((mu-a)/sigma) + p * (Phi((b-mu)/sigma) - Phi((a-mu)/sigma)))
*
*  where phi() is PDF for N(0, 1) and Phi() is CDF for N(0, 1)
*
*/


double pdf_tnorm(double x, double mu, double sigma, double a, double b) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(a) || ISNAN(b))
    return NA_REAL;
  if (sigma <= 0.0 || b <= a) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (a == R_NegInf && b == R_PosInf)
    return R::dnorm(x, mu, sigma, false);
  
  double Phi_a, Phi_b;
  if (x > a && x < b) {
    Phi_a = Phi((a-mu)/sigma);
    Phi_b = Phi((b-mu)/sigma);
    return exp(-pow(x-mu, 2.0) / (2.0*pow(sigma, 2.0))) / (SQRT_2_PI*sigma * (Phi_b - Phi_a));
  } else {
    return 0.0;
  }
}

double cdf_tnorm(double x, double mu, double sigma, double a, double b) {
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(a) || ISNAN(b))
    return NA_REAL;
  if (sigma <= 0.0 || b <= a) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (a == R_NegInf && b == R_PosInf)
    return R::pnorm(x, mu, sigma, true, false);
  
  double Phi_x, Phi_a, Phi_b;
  if (x > a && x < b) {
    Phi_x = Phi((x-mu)/sigma);
    Phi_a = Phi((a-mu)/sigma);
    Phi_b = Phi((b-mu)/sigma);
    return (Phi_x - Phi_a) / (Phi_b - Phi_a);
  } else if (x >= b) {
    return 1.0;
  } else {
    return 0.0;
  }
}

double invcdf_tnorm(double p, double mu, double sigma, double a, double b) {
  if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma) || ISNAN(a) || ISNAN(b))
    return NA_REAL;
  if (sigma <= 0.0 || b <= a || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (a == R_NegInf && b == R_PosInf)
    return R::qnorm(p, mu, sigma, true, false);
  
  double Phi_a, Phi_b;
  Phi_a = Phi((a-mu)/sigma);
  Phi_b = Phi((b-mu)/sigma);
  return InvPhi(Phi_a + p * (Phi_b - Phi_a)) * sigma + mu;
}

double rng_tnorm(double mu, double sigma, double a, double b) {
  if (ISNAN(mu) || ISNAN(sigma) || ISNAN(a) || ISNAN(b))
    return NA_REAL;
  if (sigma <= 0.0 || b <= a) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  // non-truncated normal
  if (a == R_NegInf && b == R_PosInf)
    return R::rnorm(mu, sigma);

  double r, u, za, zb, aa, za_sq, zb_sq;
  bool stop = false;

  za = (a-mu)/sigma;
  zb = (b-mu)/sigma;
  za_sq = pow(za, 2.0);
  zb_sq = pow(zb, 2.0);
  
  if (abs(za) <= 1e-16 && zb == R_PosInf) {
    r = R::norm_rand();
    if (r < 0.0)
      r = -r;
  } else if (za == R_PosInf && abs(zb) <= 1e-16) {
    r = R::norm_rand();
    if (r > 0.0)
      r = -r;
  } else if ((za < 0.0 && zb == R_PosInf) || (za == R_NegInf && zb > 0.0) ||
      (za != R_PosInf && zb != R_PosInf &&
       za < 0.0 && zb > 0.0 && zb-za > SQRT_2_PI)) {
    while (!stop) {
      r = R::norm_rand();
      if (r >= za && r <= zb)
        stop = true;
    }
  } else if (za >= 0.0 && (zb > za + 2.0*sqrt(M_E) / (za + sqrt(za_sq + 4.0))
                      * exp((za*2.0 - za*sqrt(za_sq + 4.0)) / 4.0))) {
    aa = (za + sqrt(za_sq + 4.0)) / 2.0;
    while (!stop) {
      r = R::exp_rand() / aa + za;
      u = rng_unif();
      if ((u <= exp(-pow(r-aa, 2.0) / 2.0)) && (r <= zb))
        stop = true;
    }
  } else if (zb <= 0.0 && (-za > -zb + 2.0*sqrt(M_E) / (-zb + sqrt(zb_sq + 4.0))
                          * exp((zb*2.0 + zb*sqrt(zb_sq + 4.0)) / 4.0))) {
    aa = (-zb + sqrt(zb_sq + 4.0)) / 2.0;
    while (!stop) {
      r = R::exp_rand() / aa - zb;
      u = rng_unif();
      if ((u <= exp(-pow(r-aa, 2.0) / 2.0)) && (r >= za)) {
        r = -r;
        stop = true;
      }
    }
  } else {
    if (0.0 < za) {
      while (!stop) {
        r = R::runif(za, zb);
        u = rng_unif();
        stop = (u <= exp((za_sq - pow(r, 2.0))/2.0));
      }
    } else if (zb < 0.0) {
      while (!stop) {
        r = R::runif(za, zb);
        u = rng_unif();
        stop = (u <= exp((zb_sq - pow(r, 2.0))/2.0));
      }
    } else {
      while (!stop) {
        r = R::runif(za, zb);
        u = rng_unif();
        stop = (u <= exp(-pow(r, 2.0)/2.0));
      }
    }
  }

  return mu + sigma * r;
}


// [[Rcpp::export]]
NumericVector cpp_dtnorm(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& a,
    const NumericVector& b,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_tnorm(x[i % dims[0]], mu[i % dims[1]], sigma[i % dims[2]],
                     a[i % dims[3]], b[i % dims[4]]);

  if (log_prob)
    p = Rcpp::log(p);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_ptnorm(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& a,
    const NumericVector& b,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_tnorm(x[i % dims[0]], mu[i % dims[1]], sigma[i % dims[2]],
                     a[i % dims[3]], b[i % dims[4]]);

  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qtnorm(
    const NumericVector& p,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& a,
    const NumericVector& b,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;

  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_tnorm(pp[i % dims[0]], mu[i % dims[1]], sigma[i % dims[2]],
                        a[i % dims[3]], b[i % dims[4]]);

  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rtnorm(
    const int& n,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& a,
    const NumericVector& b
  ) {

  std::vector<int> dims;
  dims.push_back(mu.length());
  dims.push_back(sigma.length());
  dims.push_back(a.length());
  dims.push_back(b.length());
  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_tnorm(mu[i % dims[0]], sigma[i % dims[1]],
                     a[i % dims[2]], b[i % dims[3]]);

  return x;
}

