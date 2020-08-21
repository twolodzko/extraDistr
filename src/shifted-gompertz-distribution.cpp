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


/*
*  Gompertz distribution
*
*  Values:
*  x >= 0
*
*  Parameters:
*  b > 0
*  eta > 0
*
*  f(x)    = b*exp(-b*x) * exp(-eta*exp(-b*x)) * (1 + eta*(1 - exp(-b*x)))
*  F(x)    = (1-exp(-b*x)) * exp(-eta*exp(-b*x))
*
* References:
*
* Bemmaor, A.C. (1994).
* Modeling the Diffusion of New Durable Goods: Word-of-Mouth Effect Versus Consumer Heterogeneity.
* [In:] G. Laurent, G.L. Lilien & B. Pras. Research Traditions in Marketing.
* Boston: Kluwer Academic Publishers. pp. 201-223.
* 
* Jimenez, T.F., Jodra, P. (2009).
* A Note on the Moments and Computer Generation of the Shifted Gompertz Distribution.
* Communications in Statistics - Theory and Methods, 38(1), 78-89.
* 
* Jimenez T.F. (2014).
* Estimation of the Parameters of the Shifted Gompertz Distribution,
* Using Least Squares, Maximum Likelihood and Moments Methods.
* Journal of Computational and Applied Mathematics, 255(1), 867-877.
*
*/


inline double logpdf_sgomp(double x, double b, double eta,
                           bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(b) || ISNAN(eta))
    return x+b+eta;
#endif
  if (b <= 0.0 || eta <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0 || !R_FINITE(x))
    return R_NegInf;
  double ebx = exp(-b*x);
  // b*ebx * exp(-eta*ebx) * (1+eta*(1-ebx));
  return log(b) + log(ebx) - eta*ebx + log1p(eta*(1-ebx));
}

inline double cdf_sgomp(double x, double b, double eta,
                           bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(b) || ISNAN(eta))
    return x+b+eta;
#endif
  if (b <= 0.0 || eta <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  if (x == R_PosInf)
    return 1.0;
  double ebx = exp(-b*x);
  // (1-ebx) * exp(-eta*ebx)
  return exp(log1p(-ebx) - eta*ebx);
}

inline double rng_sgomp(double b, double eta, bool& throw_warning) {
  if (ISNAN(b) || ISNAN(eta) || b <= 0.0 || eta <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  double u, v, rg, re;
  u = R::exp_rand(); // -log(rng_unif())
  v = R::exp_rand(); // -log(rng_unif())
  rg = -log(u/eta) / b;
  re = v / b;
  return (rg>re) ? rg : re;
}


// [[Rcpp::export]]
NumericVector cpp_dsgomp(
    const NumericVector& x,
    const NumericVector& b,
    const NumericVector& eta,
    bool log_prob = false
  ) {
  
  if (std::min({x.length(), b.length(), eta.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    b.length(),
    eta.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_sgomp(GETV(x, i), GETV(b, i),
                        GETV(eta, i), throw_warning);
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_psgomp(
    const NumericVector& x,
    const NumericVector& b,
    const NumericVector& eta,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), b.length(), eta.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    b.length(),
    eta.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_sgomp(GETV(x, i), GETV(b, i),
                     GETV(eta, i), throw_warning);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rsgomp(
    const int& n,
    const NumericVector& b,
    const NumericVector& eta
  ) {
  
  if (std::min({b.length(), eta.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_sgomp(GETV(b, i), GETV(eta, i),
                     throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

