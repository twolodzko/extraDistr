#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


/*
*  Bivariate Normal distribution
*
*  Values:
*  x
*
*  Parameters:
*  mu1, mu2
*  sigma1, sigma2 > 0
*
*  z1 = (x1 - mu1)/sigma1
*  z2 = (x2 - mu2)/sigma2
*
*  f(x)    = 1/(2*pi*sqrt(1-rho^2)*sigma1*sigma2) * exp(-(1/(2*(1-rho^2)*(z1^2 - 2*rho*z1*z2 + z2^2))))
*
*/


// [[Rcpp::export]]
NumericVector cpp_dbnorm(NumericMatrix x,
                         double mu1 = 0, double mu2 = 0,
                         double sigma1 = 1, double sigma2 = 1,
                         double rho = 0,
                         bool log_prob = false) {

  if (x.ncol() < 2)
    throw Rcpp::exception("Matrix q has to have two columns.");
  if (x.ncol() > 2)
    Rcpp::warning("Only first two columns of matrix x are going to be used.");
  if (sigma1 <= 0 || sigma2 <= 0)
    throw Rcpp::exception("Values of sigma1 and sigma2 should be > 0.");
  if (rho < -1 || rho > 1)
    throw Rcpp::exception("Value of rho should fit the [-1, 1] interval.");

  if (rho == 1) rho = 1 - 1e-16;
  if (rho == -1) rho = -1 + 1e-16;

  int n  = x.nrow();
  NumericVector p(n);
  double z1, z2;
  double c1 = 1/(2*M_PI*sqrt(1-pow(rho, 2))*sigma1*sigma2);
  double c2 = -1/(2*(1-pow(rho, 2)));

  for (int i = 0; i < n; i++) {
    z1 = (x(i, 0) - mu1)/sigma1;
    z2 = (x(i, 1) - mu2)/sigma2;
    p[i] = c1 * exp(c2 * (pow(z1, 2) - 2*rho*z1*z2 + pow(z2, 2)));
  }

  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pbnorm(NumericMatrix x,
                         double mu1 = 0, double mu2 = 0,
                         double sigma1 = 1, double sigma2 = 1,
                         double rho = 0,
                         bool lower_tail = true, bool log_prob = false,
                         int nsim = 10000) {

  if (x.ncol() < 2)
    throw Rcpp::exception("Matrix q has to have two columns.");
  if (x.ncol() > 2)
    Rcpp::warning("Only first two columns of matrix q are going to be used.");
  if (sigma1 <= 0 || sigma2 <= 0)
    throw Rcpp::exception("Values of sigma1 and sigma2 should be > 0.");
  if (rho < -1 || rho > 1)
    throw Rcpp::exception("Value or rho should fit the [-1, 1] interval.");

  if (rho == 1) rho = 1 - 1e-16;
  if (rho == -1) rho = -1 + 1e-16;

  int n = x.nrow();
  NumericVector p(n, 0.0);
  double u, v, y1, y2;

  for (int i = 0; i < n*nsim; i++) {
    u = R::rnorm(0, 1);
    v = R::rnorm(0, 1);
    y1 = mu1 + sigma1 * u;
    y2 = mu2 + sigma2 * (rho*u + sqrt(1-pow(rho, 2))*v);
    p[i % n] += (y1 <= x(i % n, 0) && y2 <= x(i % n, 1)) ? 1.0/nsim : 0.0;
  }

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericMatrix cpp_rbnorm(int n,
                         double mu1 = 0, double mu2 = 0,
                         double sigma1 = 1, double sigma2 = 1,
                         double rho = 0) {

  if (sigma1 <= 0 || sigma2 <= 0)
    throw Rcpp::exception("Values of sigma1 and sigma2 should be > 0.");
  if (rho < -1 || rho > 1)
    throw Rcpp::exception("Value or rho should fit the [-1, 1] interval.");

  if (rho == 1) rho = 1 - 1e-16;
  if (rho == -1) rho = -1 + 1e-16;

  NumericMatrix x(n, 2);
  double u, v;

  if (rho != 0) {
    for (int i = 0; i < n; i++) {
      u = R::rnorm(0, 1);
      v = R::rnorm(0, 1);
      x(i, 0) = mu1 + sigma1 * u;
      x(i, 1) = mu2 + sigma2 * (rho*u + sqrt(1-pow(rho, 2))*v);
    }
  } else {
    for (int i = 0; i < n; i++) {
      x(i, 0) = R::rnorm(mu1, sigma1);
      x(i, 1) = R::rnorm(mu2, sigma2);
    }
  }

  return x;
}

