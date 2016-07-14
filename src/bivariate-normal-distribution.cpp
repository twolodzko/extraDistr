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
  
  if (sigma1 <= 0.0 || sigma2 <= 0.0 || rho <= -1.0 || rho >= 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  
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
    bool log_prob = false
  ) {

  int nx  = x.length();
  int ny  = y.length();
  int nm1 = mu1.length();
  int nm2 = mu2.length();
  int ns1 = sigma1.length();
  int ns2 = sigma2.length();
  int nr = rho.length();
  int Nmax = Rcpp::max(IntegerVector::create(nx, ny, nm1, nm2, ns1, ns2, nr));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_bnorm(x[i % nx], y[i % ny],
                     mu1[i % nm1], mu2[i % nm2],
                     sigma1[i % ns1], sigma2[i % ns2],
                     rho[i % nr]);

  if (log_prob)
    for (int i = 0; i < nx; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericMatrix cpp_rbnorm(
    const int n,
    const NumericVector& mu1,
    const NumericVector& mu2,
    const NumericVector& sigma1,
    const NumericVector& sigma2,
    const NumericVector& rho
  ) {

  int nm1 = mu1.length();
  int nm2 = mu2.length();
  int ns1 = sigma1.length();
  int ns2 = sigma2.length();
  int nr = rho.length();
  
  NumericMatrix x(n, 2);

  for (int i = 0; i < n; i++) {
    if (sigma1[i % ns1] <= 0.0 || sigma2[i % ns2] <= 0.0 ||
        rho[i % nr] < -1.0 || rho[i % nr] > 1.0) {
      Rcpp::warning("NaNs produced");
      x(i, 0) = NAN;
      x(i, 1) = NAN;
    } else if (rho[i % nr] != 0.0) {
      double u = R::norm_rand();
      double v = R::norm_rand();
      double corr = (rho[i % nr]*u + sqrt(1.0 - pow(rho[i % nr], 2.0))*v);
      x(i, 0) = mu1[i % nm1] + sigma1[i % ns1] * u;
      x(i, 1) = mu2[i % nm2] + sigma2[i % ns2] * corr;
    } else {
      x(i, 0) = R::rnorm(mu1[i % nm1], sigma1[i % ns1]);
      x(i, 1) = R::rnorm(mu2[i % nm2], sigma2[i % ns2]);
    }
  }

  return x;
}

