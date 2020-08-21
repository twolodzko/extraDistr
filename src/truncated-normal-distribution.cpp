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


double pdf_tnorm(double x, double mu, double sigma,
                 double a, double b, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(a) || ISNAN(b))
    return x+mu+sigma+a+b;
#endif
  if (sigma <= 0.0 || b <= a) {
    throw_warning = true;
    return NAN;
  }
  
  if (a == R_NegInf && b == R_PosInf)
    return R::dnorm(x, mu, sigma, false);
  
  double Phi_a, Phi_b;
  if (x > a && x < b) {
    Phi_a = Phi((a-mu)/sigma);
    Phi_b = Phi((b-mu)/sigma);
    return exp(-((x-mu)*(x-mu)) / (2.0*(sigma*sigma))) /
              (SQRT_2_PI*sigma * (Phi_b - Phi_a));
  } else {
    return 0.0;
  }
}

double cdf_tnorm(double x, double mu, double sigma,
                 double a, double b, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(a) || ISNAN(b))
    return x+mu+sigma+a+b;
#endif
  if (sigma <= 0.0 || b <= a) {
    throw_warning = true;
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

double invcdf_tnorm(double p, double mu, double sigma,
                    double a, double b, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma) || ISNAN(a) || ISNAN(b))
    return p+mu+sigma+a+b;
#endif
  if (sigma <= 0.0 || b <= a || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  
  if (a == R_NegInf && b == R_PosInf)
    return R::qnorm(p, mu, sigma, true, false);
  
  double Phi_a, Phi_b;
  Phi_a = Phi((a-mu)/sigma);
  Phi_b = Phi((b-mu)/sigma);
  return InvPhi(Phi_a + p * (Phi_b - Phi_a)) * sigma + mu;
}

double rng_tnorm(double mu, double sigma, double a,
                 double b, bool& throw_warning) {
  if (ISNAN(mu) || ISNAN(sigma) || ISNAN(a) || ISNAN(b) ||
      sigma <= 0.0 || b <= a) {
    throw_warning = true;
    return NA_REAL;
  }
  
  // non-truncated normal
  if (a == R_NegInf && b == R_PosInf)
    return R::rnorm(mu, sigma);

  double r, u, za, zb, aa, za_sq, zb_sq;
  bool stop = false;

  za = (a-mu)/sigma;
  zb = (b-mu)/sigma;
  za_sq = za * za;
  zb_sq = zb * zb;
  
  if (abs(za) <= 1e-16 && zb == R_PosInf) {
    r = R::norm_rand();
    if (r < 0.0)
      r = -r;
  } else if (za == R_PosInf && abs(zb) <= 1e-16) {
    r = R::norm_rand();
    if (r > 0.0)
      r = -r;
  } else if ((za < 0.0 && zb == R_PosInf) ||
      (za == R_NegInf && zb > 0.0) ||
      (za != R_PosInf && zb != R_PosInf &&
       za < 0.0 && zb > 0.0 && zb-za > SQRT_2_PI)) {
    do {
      r = R::norm_rand();
      if (r >= za && r <= zb)
        stop = true;
    } while (!stop);
  } else if (za >= 0.0 && (zb > za + 2.0*sqrt(M_E) / (za + sqrt(za_sq + 4.0))
                      * exp((za*2.0 - za*sqrt(za_sq + 4.0)) / 4.0))) {
    aa = (za + sqrt(za_sq + 4.0)) / 2.0;
    do {
      r = R::exp_rand() / aa + za;
      u = rng_unif();
      if ((u <= exp(-((r-aa)*(r-aa)) / 2.0)) && (r <= zb))
        stop = true;
    } while (!stop);
  } else if (zb <= 0.0 && (-za > -zb + 2.0*sqrt(M_E) / (-zb + sqrt(zb_sq + 4.0))
                          * exp((zb*2.0 + zb*sqrt(zb_sq + 4.0)) / 4.0))) {
    aa = (-zb + sqrt(zb_sq + 4.0)) / 2.0;
    do {
      r = R::exp_rand() / aa - zb;
      u = rng_unif();
      if ((u <= exp(-((r-aa)*(r-aa)) / 2.0)) && (r <= -za)) {
        r = -r;
        stop = true;
      }
    } while (!stop);
  } else {
    if (0.0 < za) {
      do {
        r = R::runif(za, zb);
        u = rng_unif();
        stop = (u <= exp((za_sq - r*r)/2.0));
      } while (!stop);
    } else if (zb < 0.0) {
      do {
        r = R::runif(za, zb);
        u = rng_unif();
        stop = (u <= exp((zb_sq - r*r)/2.0));
      } while (!stop);
    } else {
      do {
        r = R::runif(za, zb);
        u = rng_unif();
        stop = (u <= exp(-(r*r)/2.0));
      } while (!stop);
    }
  }

  return mu + sigma * r;
}


// [[Rcpp::export]]
NumericVector cpp_dtnorm(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& lower,
    const NumericVector& upper,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), mu.length(), sigma.length(),
                lower.length(), upper.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    mu.length(),
    sigma.length(),
    lower.length(),
    upper.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_tnorm(GETV(x, i), GETV(mu, i),
                     GETV(sigma, i), GETV(lower, i),
                     GETV(upper, i), throw_warning);

  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_ptnorm(
    const NumericVector& x,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& lower,
    const NumericVector& upper,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), mu.length(), sigma.length(),
                lower.length(), upper.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    mu.length(),
    sigma.length(),
    lower.length(),
    upper.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_tnorm(GETV(x, i), GETV(mu, i),
                     GETV(sigma, i), GETV(lower, i),
                     GETV(upper, i), throw_warning);

  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qtnorm(
    const NumericVector& p,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& lower,
    const NumericVector& upper,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), mu.length(), sigma.length(),
                lower.length(), upper.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    p.length(),
    mu.length(),
    sigma.length(),
    lower.length(),
    upper.length()
  });
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;

  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_tnorm(GETV(pp, i), GETV(mu, i),
                        GETV(sigma, i), GETV(lower, i),
                        GETV(upper, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rtnorm(
    const int& n,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& lower,
    const NumericVector& upper
  ) {

  if (std::min({mu.length(), sigma.length(),
                lower.length(), upper.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }

  NumericVector x(n);
  
  bool throw_warning = false;

  for (int i = 0; i < n; i++)
    x[i] = rng_tnorm(GETV(mu, i), GETV(sigma, i),
                     GETV(lower, i), GETV(upper, i),
                     throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

