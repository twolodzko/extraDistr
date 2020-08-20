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


inline double logpdf_huber(double x, double mu, double sigma,
                           double c, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(c))
    return x+mu+sigma+c;
#endif
  if (sigma <= 0.0 || c <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  
  double z, A, rho;
  z = abs((x - mu)/sigma);
  // A = 2.0*SQRT_2_PI * (Phi(c) + phi(c)/c - 0.5);
  A = LOG_2F + log(SQRT_2_PI) + log(Phi(c) + phi(c)/c - 0.5);

  if (z <= c) {
    rho = (z*z)/2.0;
  } else {
    rho = c*z - (c*c)/2.0;
  }

  // exp(-rho)/A/sigma;
  return -rho - A - log(sigma);
}

inline double cdf_huber(double x, double mu, double sigma,
                        double c, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma) || ISNAN(c))
    return x+mu+sigma+c;
#endif
  if (sigma <= 0.0 || c <= 0.0) {
    throw_warning = true;
    return NAN;
  }

  double A, z, az, p;
  A = 2.0*(phi(c)/c - Phi(-c) + 0.5);
  z = (x - mu)/sigma;
  az = -abs(z);
  
  if (az <= -c) 
    p = exp((c*c)/2.0)/c * exp(c*az) / SQRT_2_PI/A;
  else
    p = (phi(c)/c + Phi(az) - Phi(-c))/A;
  
  if (z <= 0.0)
    return p;
  else
    return 1.0 - p;
}

inline double invcdf_huber(double p, double mu, double sigma,
                           double c, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma) || ISNAN(c))
    return p+mu+sigma+c;
#endif
  if (sigma <= 0.0 || c <= 0.0 || !VALID_PROB(p)) {
    throw_warning = true;
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

inline double rng_huber(double mu, double sigma, double c,
                        bool& throw_warning) {
  if (ISNAN(mu) || ISNAN(sigma) || ISNAN(c) ||
      sigma <= 0.0 || c <= 0.0) {
    throw_warning = true;
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
  
  if (std::min({x.length(), mu.length(),
                sigma.length(), epsilon.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    mu.length(),
    sigma.length(),
    epsilon.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_huber(GETV(x, i), GETV(mu, i),
                        GETV(sigma, i), GETV(epsilon, i),
                        throw_warning);
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
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
  
  if (std::min({x.length(), mu.length(),
                sigma.length(), epsilon.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    mu.length(),
    sigma.length(),
    epsilon.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_huber(GETV(x, i), GETV(mu, i),
                     GETV(sigma, i), GETV(epsilon, i),
                     throw_warning);
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
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
  
  if (std::min({p.length(), mu.length(),
                sigma.length(), epsilon.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    p.length(),
    mu.length(),
    sigma.length(),
    epsilon.length()
  });
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_huber(GETV(pp, i), GETV(mu, i),
                        GETV(sigma, i), GETV(epsilon, i),
                        throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rhuber(
    const int& n,
    const NumericVector& mu,
    const NumericVector& sigma,
    const NumericVector& epsilon
  ) {
  
  if (std::min({mu.length(), sigma.length(), epsilon.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_huber(GETV(mu, i), GETV(sigma, i),
                     GETV(epsilon, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

