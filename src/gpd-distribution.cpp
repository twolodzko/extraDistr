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
*  Generalized Pareto distribution
*
*  Values:
*  x
*
*  Parameters:
*  mu
*  sigma > 0
*  xi
*
*  z = (x-mu)/sigma
*  where 1+xi*z > 0
*
*  f(x)    = { (1+xi*z)^{-(xi+1)/xi}/sigma       if xi != 0
*            { exp(-z)/sigma                     otherwise
*  F(x)    = { 1-(1+xi*z)^{-1/xi}                if xi != 0
*            { 1-exp(-z)                         otherwise
*  F^-1(p) = { mu + sigma * ((1-p)^{-xi}-1)/xi   if xi != 0
*            { mu - sigma * log(1-p)             otherwise
*
*/

double pdf_gpd(double x, double mu, double sigma, double xi) {
  if (sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double z = (x-mu)/sigma;
  if (xi != 0.0) {
    if (x >= mu)
      return pow(1.0+xi*z, -(xi+1.0)/xi)/sigma;
    else
      return 0.0;
  } else {
    if (x >= mu && x <= (mu - sigma/xi))
      return exp(-z)/sigma;
    else
      return 0.0;
  }
}

double cdf_gpd(double x, double mu, double sigma, double xi) {
  if (sigma <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double z = (x-mu)/sigma;
  if (xi != 0.0) {
    if (x >= mu)
      return 1.0 - pow(1.0+xi*z, -1.0/xi);
    else
      return 0.0;
  } else {
    if (x >= mu && x <= (mu - sigma/xi))
      return 1.0 - exp(-z);
    else
      return 0.0;
  }
}

double invcdf_gpd(double p, double mu, double sigma, double xi) {
  if (sigma <= 0.0 || p < 0.0 || p > 1.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (xi != 0.0)
    return mu + sigma * (pow(1.0-p, -xi)-1.0)/xi;
  else
    return mu - sigma * log(1.0-p);
}


// [[Rcpp::export]]
NumericVector cpp_dgpd(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& xi,
    bool log_prob = false
  ) {

  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int nx = xi.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, nx));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_gpd(x[i % n], mu[i % nm], sigma[i % ns], xi[i % nx]);

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pgpd(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& xi,
    bool lower_tail = true, bool log_prob = false
  ) {

  int n  = x.length();
  int nm = mu.length();
  int ns = sigma.length();
  int nx = xi.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, nx));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_gpd(x[i % n], mu[i % nm], sigma[i % ns], xi[i % nx]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1.0 - p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qgpd(
    const NumericVector& p,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& xi,
    bool lower_tail = true, bool log_prob = false
  ) {

  int n  = p.length();
  int nm = mu.length();
  int ns = sigma.length();
  int nx = xi.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, ns, nx));
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);

  if (log_prob)
    for (int i = 0; i < n; i++)
      pp[i] = exp(p[i]);

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      pp[i] = 1.0 - pp[i];

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_gpd(pp[i % n], mu[i % nm], sigma[i % ns], xi[i % nx]);

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rgpd(
    const int n,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& xi
  ) {

  double u;
  int nm = mu.length();
  int ns = sigma.length();
  int nx = xi.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    u = R::runif(0.0, 1.0);
    x[i] = invcdf_gpd(u, mu[i % nm], sigma[i % ns], xi[i % nx]);
  }

  return x;
}

