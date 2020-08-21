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
Joiner, B.L., & Rosenblatt, J.R. (1971).
Some properties of the range in samples from Tukey's symmetric lambda distributions.
Journal of the American Statistical Association, 66(334), 394-399.

Hastings Jr, C., Mosteller, F., Tukey, J.W., & Winsor, C.P. (1947).
Low moments for small samples: a comparative study of order statistics.
The Annals of Mathematical Statistics, 413-426.
*/


inline double invcdf_tlambda(double p, double lambda,
                             bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(p) || ISNAN(lambda))
    return p+lambda;
#endif
  if (!VALID_PROB(p)) {
    throw_warning = true;
    return NAN;
  }
  if (lambda == 0.0)
    return log(p) - log(1.0 - p);
  return (pow(p, lambda) - pow(1.0 - p, lambda))/lambda;
}

inline double rng_tlambda(double lambda, bool& throw_warning) {
  if (ISNAN(lambda)) {
    throw_warning = true;
    return NA_REAL;
  }
  double u = rng_unif();
  if (lambda == 0.0)
    return log(u) - log(1.0 - u);
  return (pow(u, lambda) - pow(1.0 - u, lambda))/lambda;
}


// [[Rcpp::export]]
NumericVector cpp_qtlambda(
    const NumericVector& p,
    const NumericVector& lambda,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  if (std::min({p.length(), lambda.length()}) < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max({
    p.length(),
    lambda.length()
  });
  NumericVector q(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_tlambda(GETV(pp, i), GETV(lambda, i),
                          throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rtlambda(
    const int& n,
    const NumericVector& lambda
  ) {
  
  if (lambda.length() < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }
  
  NumericVector x(n);
  
  bool throw_warning = false;
    
  for (int i = 0; i < n; i++)
    x[i] = rng_tlambda(GETV(lambda, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

