#include <Rcpp.h>
#include "shared.h"
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
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


inline double pdf_bnorm(double x, double y, double mu1, double mu2,
                        double sigma1, double sigma2, double rho,
                        bool& throw_warning) {
  
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(y) || ISNAN(mu1) || ISNAN(mu2) ||
      ISNAN(sigma1) || ISNAN(sigma2) || ISNAN(rho))
    return x+y+mu1+mu2+sigma1+sigma2+rho;
#endif
  
  if (sigma1 <= 0.0 || sigma2 <= 0.0 || rho <= -1.0 || rho >= 1.0) {
    throw_warning = true;
    return NAN;
  }
  
  if (!R_FINITE(x) || !R_FINITE(y))
    return 0.0;
  
  double z1 = (x - mu1)/sigma1;
  double z2 = (y - mu2)/sigma2;
  
  double c1 = 1.0/(2.0*M_PI*sqrt(1.0 - (rho*rho))*sigma1*sigma2);
  double c2 = -1.0/(2.0*(1.0 - (rho*rho)));
  
  return c1 * exp(c2 * ((z1*z1) - 2.0*rho*z1*z2 + (z2*z2)));
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
  
  if (std::min({x.length(), y.length(),
                mu1.length(), mu2.length(),
                sigma1.length(), sigma2.length(),
                rho.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    y.length(),
    mu1.length(),
    mu2.length(),
    sigma1.length(),
    sigma2.length(),
    rho.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  if (x.length() != y.length())
    Rcpp::stop("lengths of x and y differ");

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_bnorm(GETV(x, i), GETV(y, i), GETV(mu1, i),
                     GETV(mu2, i), GETV(sigma1, i),
                     GETV(sigma2, i), GETV(rho, i),
                     throw_warning);

  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

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
  
  if (std::min({mu1.length(), mu2.length(),
                sigma1.length(), sigma2.length(),
                rho.length()}) < 1) {
    Rcpp::warning("NAs produced");
    NumericMatrix out(n, 2);
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }

  NumericMatrix x(n, 2);
  double u, v, corr;
  
  bool throw_warning = false;

  for (int i = 0; i < n; i++) {
    if (ISNAN(GETV(mu1, i)) || ISNAN(GETV(mu2, i)) ||
        ISNAN(GETV(sigma1, i)) || ISNAN(GETV(sigma2, i)) ||
        ISNAN(GETV(rho, i)) || GETV(sigma1, i) <= 0.0 ||
        GETV(sigma2, i) <= 0.0 || GETV(rho, i) < -1.0 ||
        GETV(rho, i) > 1.0) {
      throw_warning = true;
      x(i, 0) = NA_REAL;
      x(i, 1) = NA_REAL;
    } else if (!tol_equal(GETV(rho, i), 0.0)) {
      u = R::norm_rand();
      v = R::norm_rand();
      corr = (GETV(rho, i)*u + sqrt(1.0 - pow(GETV(rho, i), 2.0))*v);
      x(i, 0) = GETV(mu1, i) + GETV(sigma1, i) * u;
      x(i, 1) = GETV(mu2, i) + GETV(sigma2, i) * corr;
    } else {
      x(i, 0) = R::rnorm(GETV(mu1, i), GETV(sigma1, i));
      x(i, 1) = R::rnorm(GETV(mu2, i), GETV(sigma2, i));
    }
  }
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

