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
 * Birnbaum-Saunders (Fatigue Life) Distribution
 * 
 * Support:
 * x > mu
 * 
 * Parameters:
 * mu
 * alpha > 0
 * beta > 0
 * 
 * 
 */

inline double logpdf_fatigue(double x, double alpha, double beta,
                             double mu, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(alpha) || ISNAN(beta) || ISNAN(mu))
    return x+alpha+beta+mu;
#endif
  if (alpha <= 0.0 || beta <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x <= mu || !R_FINITE(x))
    return R_NegInf;
  double z, zb, bz;
  z = x-mu;
  zb = sqrt(z/beta);
  bz = sqrt(beta/z);
  // (zb+bz)/(2.0*alpha*z) * phi((zb-bz)/alpha)
  return log(zb+bz) - LOG_2F - log(alpha) - log(z) + lphi((zb-bz)/alpha);
}

inline double cdf_fatigue(double x, double alpha, double beta,
                          double mu, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(alpha) || ISNAN(beta) || ISNAN(mu))
    return x+alpha+beta+mu;
#endif
  if (alpha <= 0.0 || beta <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x <= mu)
    return 0.0;
  double z, zb, bz;
  z = x-mu;
  zb = sqrt(z/beta);
  bz = sqrt(beta/z);
  return Phi((zb-bz)/alpha);
}

inline double invcdf_fatigue(double p, double alpha, double beta,
                             double mu, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(alpha) || ISNAN(beta) || ISNAN(mu))
    return p+alpha+beta+mu;
#endif
  if (alpha <= 0.0 || beta <= 0.0 || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  if (p == 0.0)
    return mu;
  double Zp = InvPhi(p);
  return pow(alpha/2.0*Zp + sqrt(pow(alpha/2.0*Zp, 2.0) + 1.0), 2.0) * beta + mu;
}

inline double rng_fatigue(double alpha, double beta,
                          double mu, bool& throw_warning) {
  if (ISNAN(alpha) || ISNAN(beta) || ISNAN(mu) || alpha <= 0.0 || beta <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  double z = R::norm_rand();
  return pow(alpha/2.0*z + sqrt(pow(alpha/2.0*z, 2.0) + 1.0), 2.0) * beta + mu;
}


// [[Rcpp::export]]
NumericVector cpp_dfatigue(
    const NumericVector& x,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& mu,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), alpha.length(),
                beta.length(), mu.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    alpha.length(),
    beta.length(),
    mu.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_fatigue(GETV(x, i), GETV(alpha, i),
                          GETV(beta, i), GETV(mu, i),
                          throw_warning);
  
  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pfatigue(
    const NumericVector& x,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& mu,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), alpha.length(),
                beta.length(), mu.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    alpha.length(),
    beta.length(),
    mu.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_fatigue(GETV(x, i), GETV(alpha, i),
                       GETV(beta, i), GETV(mu, i),
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
NumericVector cpp_qfatigue(
    const NumericVector& p,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& mu,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), alpha.length(),
                beta.length(), mu.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    p.length(),
    alpha.length(),
    beta.length(),
    mu.length()
  });
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_fatigue(GETV(pp, i), GETV(alpha, i),
                          GETV(beta, i), GETV(mu, i),
                          throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rfatigue(
    const int& n,
    const NumericVector& alpha,
    const NumericVector& beta,
    const NumericVector& mu
  ) {
  
  if (std::min({alpha.length(), beta.length(), mu.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_fatigue(GETV(alpha, i), GETV(beta, i),
                       GETV(mu, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

