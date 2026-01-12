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
 *  Nakagami distribution
 *
 *  Values:
 *  x >= 0
 *
 *  Parameters:
 *  m > 0     (shape)
 *  w > 0     (scale, omega)
 *
 *  f(x)    = 2*m^m / (gamma(m)*w^m) * x^(2m-1) * exp(-m/w * x^2)
 *  F(x)    = pgamma(x^2, shape = m, rate = m/w)
 *  F^-1(p) = sqrt(qgamma(p, shape = m, rate = m/w))
 *
 *  If Y ~ Gamma(shape = m, rate = m/w), then X = sqrt(Y) ~ Nakagami(m, w)
 *
 */


inline double logpdf_nakagami(double x, double m, double w,
                               bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(m) || ISNAN(w))
    return x+m+w;
#endif
  if (m <= 0.0 || w <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x <= 0.0 || !R_FINITE(x))
    return R_NegInf;
  // f(x) = 2 * x * dgamma(x^2, shape = m, scale = w/m)
  // log(f(x)) = log(2) + log(x) + log(dgamma(x^2, shape = m, scale = w/m))
  double scale = w / m;
  return LOG_2F + log(x) + R::dgamma(x*x, m, scale, true);
}

inline double cdf_nakagami(double x, double m, double w,
                            bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(m) || ISNAN(w))
    return x+m+w;
#endif
  if (m <= 0.0 || w <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (x <= 0.0)
    return 0.0;
  if (!R_FINITE(x))
    return 1.0;
  double scale = w / m;
  return R::pgamma(x*x, m, scale, true, false);
}

inline double invcdf_nakagami(double p, double m, double w,
                               bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(m) || ISNAN(w))
    return p+m+w;
#endif
  if (!VALID_PROB(p) || m <= 0.0 || w <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  double scale = w / m;
  return sqrt(R::qgamma(p, m, scale, true, false));
}

inline double rng_nakagami(double m, double w, bool& throw_warning) {
  if (ISNAN(m) || ISNAN(w) || m <= 0.0 || w <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  double scale = w / m;
  return sqrt(R::rgamma(m, scale));
}


// [[Rcpp::export]]
NumericVector cpp_dnaka(
    const NumericVector& x,
    const NumericVector& m,
    const NumericVector& w,
    const bool& log_prob = false
  ) {

  if (std::min({x.length(), m.length(), w.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    m.length(),
    w.length()
  });
  NumericVector p(Nmax);

  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = logpdf_nakagami(GETV(x, i), GETV(m, i), GETV(w, i),
                           throw_warning);

  if (!log_prob)
    p = Rcpp::exp(p);

  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pnaka(
    const NumericVector& x,
    const NumericVector& m,
    const NumericVector& w,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  if (std::min({x.length(), m.length(), w.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    m.length(),
    w.length()
  });
  NumericVector p(Nmax);

  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_nakagami(GETV(x, i), GETV(m, i), GETV(w, i),
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
NumericVector cpp_qnaka(
    const NumericVector& p,
    const NumericVector& m,
    const NumericVector& w,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  if (std::min({p.length(), m.length(), w.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    p.length(),
    m.length(),
    w.length()
  });
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);

  bool throw_warning = false;

  if (log_prob)
    pp = Rcpp::exp(pp);

  if (!lower_tail)
    pp = 1.0 - pp;

  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_nakagami(GETV(pp, i), GETV(m, i), GETV(w, i),
                           throw_warning);

  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rnaka(
    const int& n,
    const NumericVector& m,
    const NumericVector& w
  ) {

  if (std::min({m.length(), w.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }

  NumericVector x(n);

  bool throw_warning = false;

  for (int i = 0; i < n; i++)
    x[i] = rng_nakagami(GETV(m, i), GETV(w, i), throw_warning);

  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}
