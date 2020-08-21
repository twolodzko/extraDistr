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
*  Bernoulli distribution
*
*  Values:
*  x
*
*  Parameters:
*  0 <= p <= 1
*
*/

inline double pdf_bernoulli(double x, double prob,
                            bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(prob))
    return x+prob;
#endif
  if (!VALID_PROB(prob)) {
    throw_warning = true;
    return NAN;
  }
  if (x == 1.0)
    return prob;
  if (x == 0.0)
    return 1.0 - prob;
  
  char msg[55];
  std::snprintf(msg, sizeof(msg), "improper x = %f", x);
  Rcpp::warning(msg);
  
  return 0.0;
}

inline double cdf_bernoulli(double x, double prob,
                            bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(prob))
    return x+prob;
#endif
  if (!VALID_PROB(prob)) {
    throw_warning = true;
    return NAN;
  }
  if (x < 0.0)
    return 0.0;
  if (x < 1.0)
    return 1.0 - prob;
  return 1.0;
}

inline double invcdf_bernoulli(double p, double prob,
                               bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(prob))
    return p+prob;
#endif
  if (!VALID_PROB(prob) || !VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  return (p <= (1.0 - prob)) ? 0.0 : 1.0;
}

inline double rng_bernoulli(double prob, bool& throw_warning) {
  if (ISNAN(prob) || !VALID_PROB(prob)) {
    throw_warning = true;
    return NA_REAL;
  }
  double u = rng_unif();
  return (u > prob) ? 0.0 : 1.0;
}


// [[Rcpp::export]]
NumericVector cpp_dbern(
    const NumericVector& x,
    const NumericVector& prob,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), prob.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    prob.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_bernoulli(GETV(x, i), GETV(prob, i),
                         throw_warning);
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pbern(
    const NumericVector& x,
    const NumericVector& prob,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({x.length(), prob.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    x.length(),
    prob.length()
  });
  NumericVector p(Nmax);

  bool throw_warning = false;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_bernoulli(GETV(x, i), GETV(prob, i),
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
NumericVector cpp_qbern(
    const NumericVector& p,
    const NumericVector& prob,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), prob.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    p.length(),
    prob.length()
  });
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_bernoulli(GETV(pp, i), GETV(prob, i),
                            throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rbern(
    const int& n,
    const NumericVector& prob
  ) {
  
  if (prob.length() < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
  
  for (int i = 0; i < n; i++)
    x[i] = rng_bernoulli(GETV(prob, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

