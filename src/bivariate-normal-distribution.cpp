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
*  Bivariate Normal distribution
*
*  Values:
*  x, y
*
*  Parameters:
*  mu1, mu2
*  sigma1, sigma2 > 0
*
*  z1 = (x1 - mu1)/sigma1
*  z2 = (x2 - mu2)/sigma2
*
*  f(x) = 1/(2*pi*sqrt(1-rho^2)*sigma1*sigma2) *
*         exp(-(1/(2*(1-rho^2)*(z1^2 - 2*rho*z1*z2 + z2^2))))
*
*/


double pdf_bnorm(double x, double y,
                 double mu1, double mu2,
                 double sigma1, double sigma2,
                 double rho) {
  
  if (ISNAN(x) || ISNAN(y) || ISNAN(mu1) || ISNAN(mu2) ||
      ISNAN(sigma1) || ISNAN(sigma2) || ISNAN(rho))
    return NA_REAL;
  
  if (sigma1 <= 0.0 || sigma2 <= 0.0 || rho <= -1.0 || rho >= 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
  if (!R_finite(x) || !R_finite(y))
    return 0.0;
  
  double z1 = (x - mu1)/sigma1;
  double z2 = (y - mu2)/sigma2;
  
  double c1 = 1.0/(2.0*M_PI*sqrt(1.0 - pow(rho, 2.0))*sigma1*sigma2);
  double c2 = -1.0/(2.0*(1.0 - pow(rho, 2.0)));
  
  return c1 * exp(c2 * (pow(z1, 2.0) - 2.0*rho*z1*z2 + pow(z2, 2.0)));
}


// [[Rcpp::export]]
NumericVector cpp_dbnorm(
    const NumericVector& x,
    const NumericVector& y,
    const NumericVector& mu1,
    const NumericVector& mu2,
    const NumericVector& sigma1,
    const NumericVector& sigma2,
    const NumericVector& rho,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(y.length());
  dims.push_back(mu1.length());
  dims.push_back(mu2.length());
  dims.push_back(sigma1.length());
  dims.push_back(sigma2.length());
  dims.push_back(rho.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  if (dims[0] != dims[1])
    Rcpp::stop("lengths of x and y differ");

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_bnorm(x[i % dims[0]], y[i % dims[1]],
                     mu1[i % dims[2]], mu2[i % dims[3]],
                     sigma1[i % dims[4]], sigma2[i % dims[5]],
                     rho[i % dims[6]]);

  if (log_prob)
    p = Rcpp::log(p);

  return p;
}


// [[Rcpp::export]]
NumericMatrix cpp_rbnorm(
    const int& n,
    const NumericVector& mu1,
    const NumericVector& mu2,
    const NumericVector& sigma1,
    const NumericVector& sigma2,
    const NumericVector& rho
  ) {

  std::vector<int> dims;
  dims.push_back(mu1.length());
  dims.push_back(mu2.length());
  dims.push_back(sigma1.length());
  dims.push_back(sigma2.length());
  dims.push_back(rho.length());
  NumericMatrix x(n, 2);

  for (int i = 0; i < n; i++) {
    if (ISNAN(mu1[i % dims[0]]) || ISNAN(mu2[i % dims[1]]) ||
        ISNAN(sigma1[i % dims[2]]) || ISNAN(sigma2[i % dims[3]]) ||
        ISNAN(rho[i % dims[4]]) ||
        sigma1[i % dims[2]] <= 0.0 || sigma2[i % dims[3]] <= 0.0 ||
        rho[i % dims[4]] < -1.0 || rho[i % dims[4]] > 1.0) {
      Rcpp::warning("NAs produced");
      x(i, 0) = NA_REAL;
      x(i, 1) = NA_REAL;
    } else if (!tol_equal(rho[i % dims[4]], 0.0)) {
      double u = R::norm_rand();
      double v = R::norm_rand();
      double corr = (rho[i % dims[4]]*u + sqrt(1.0 - pow(rho[i % dims[4]], 2.0))*v);
      x(i, 0) = mu1[i % dims[0]] + sigma1[i % dims[2]] * u;
      x(i, 1) = mu2[i % dims[1]] + sigma2[i % dims[3]] * corr;
    } else {
      x(i, 0) = R::rnorm(mu1[i % dims[0]], sigma1[i % dims[2]]);
      x(i, 1) = R::rnorm(mu2[i % dims[1]], sigma2[i % dims[3]]);
    }
  }

  return x;
}

