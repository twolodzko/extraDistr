#include <Rcpp.h>
#include "shared.h"

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;


/*
*  Beta-negative binomial distribution
*
*  Values:
*  x
*
*  Parameters:
*  r > 0
*  alpha > 0
*  beta > 0
*
*  f(k) = gamma(r+k)/(k! gamma(r)) * beta(alpha+r, beta+k)/beta(alpha, beta)
*
*/

inline double pmf_bnbinom(double k, double r, double alpha,
                          double beta, bool& throw_warning) {
  if (ISNAN(k) || ISNAN(r) || ISNAN(alpha) || ISNAN(beta))
    return k+r+alpha+beta;
  if (alpha <= 0.0 || beta <= 0.0 || r < 0.0 || !isInteger(r, false)) {
    throw_warning = true;
    return NAN;
  }
  if (!isInteger(k) || k < 0.0 || !R_FINITE(k))
    return 0.0;
  return (R::gammafn(r+k) / (R::gammafn(k+1.0) * R::gammafn(r))) *
          R::beta(alpha+r, beta+k) / R::beta(alpha, beta);
}

inline double logpmf_bnbinom(double k, double r, double alpha,
                             double beta, bool& throw_warning) {
  if (ISNAN(k) || ISNAN(r) || ISNAN(alpha) || ISNAN(beta))
    return k+r+alpha+beta;
  if (alpha <= 0.0 || beta <= 0.0 || r < 0.0 || !isInteger(r, false)) {
    throw_warning = true;
    return NAN;
  }
  if (!isInteger(k) || k < 0.0 || !R_FINITE(k))
    return R_NegInf;
  return (R::lgammafn(r+k) - (R::lgammafn(k+1.0) + R::lgammafn(r))) +
    R::lbeta(alpha+r, beta+k) - R::lbeta(alpha, beta);
}

inline std::vector<double> cdf_bnbinom_table(double k, double r,
                                             double alpha, double beta) {
  
  if (k < 0.0 || !R_FINITE(k) || r < 0.0 || alpha < 0.0 || beta < 0.0)
    Rcpp::stop("inadmissible values");

  long int ik = to_int(k);
  std::vector<double> p_tab(ik+1);
  double grx, xf, gr, gar, gbx, gabrx, bab;
  
  bab = R::lbeta(alpha, beta);
  gr = R::lgammafn(r);
  gar = R::lgammafn(alpha + r);
  xf = 0.0;
  
  // k < 1
  
  grx = gr;
  gbx = R::lgammafn(beta);
  gabrx = R::lgammafn(alpha + beta + r);
  p_tab[0] = exp(grx - gr + gar + gbx - gabrx - bab);
  
  if (ik < 1)
    return p_tab;
  
  // k < 2
  
  grx += log(r);
  gbx += log(beta);
  gabrx += log(alpha + beta + r);
  p_tab[1] = p_tab[0] + exp(grx - gr + gar + gbx - gabrx - bab);
  
  if (ik < 2)
    return p_tab;
  
  // k >= 2
  
  double dj;
  
  for (long int j = 2; j <= ik; j++) {
    dj = to_dbl(j);
    grx += log(r + dj - 1.0);
    gbx += log(beta + dj - 1.0);
    gabrx += log(alpha + beta + r + dj - 1.0);
    xf += log(dj);
    p_tab[j] = p_tab[j-1] +
      exp(grx - (xf + gr) + gar + gbx - gabrx - bab);
  }
  
  return p_tab;
}

inline double rng_bnbinom(double r, double alpha,
                          double beta, bool& throw_warning) {
  if (ISNAN(r) || ISNAN(alpha) || ISNAN(beta) || alpha <= 0.0 ||
      beta <= 0.0 || r < 0.0 || !isInteger(r, false)) {
    throw_warning = true;
    return NA_REAL;
  }
  double prob = R::rbeta(alpha, beta);
  return R::rnbinom(r, prob);
}


// [[Rcpp::export]]
NumericVector cpp_dbnbinom(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& log_prob = false
  ) {

  int Nmax = std::max({
    x.length(),
    size.length(),
    alpha.length(),
    beta.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = logpmf_bnbinom(GETV(x, i), GETV(size, i), GETV(alpha, i),
                          GETV(beta, i), throw_warning);

  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pbnbinom(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  int Nmax = std::max({
    x.length(),
    size.length(),
    alpha.length(),
    beta.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  std::map<std::tuple<int, int, int>, std::vector<double>> memo;
  double mx = finite_max(x);
  
  for (int i = 0; i < Nmax; i++) {
    
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    
    if (ISNAN(GETV(x, i)) || ISNAN(GETV(size, i)) ||
        ISNAN(GETV(alpha, i)) || ISNAN(GETV(beta, i))) {
      p[i] = GETV(x, i) + GETV(size, i) + GETV(alpha, i) + GETV(beta, i);
    } else if (GETV(alpha, i) <= 0.0 || GETV(beta, i) <= 0.0 ||
               GETV(size, i) < 0.0 || !isInteger(GETV(size, i), false)) {
      throw_warning = true;
      p[i] = NAN;
    } else if (GETV(x, i) < 0.0) {
      p[i] = 0.0;
    } else if (!R_FINITE(GETV(x, i))) {
      p[i] = 1.0;
    } else {
      
      std::vector<double>& tmp = memo[std::make_tuple(i % size.length(),
                                                      i % alpha.length(),
                                                      i % beta.length())];
      if (!tmp.size()) {
        tmp = cdf_bnbinom_table(mx, GETV(size, i), GETV(alpha, i), GETV(beta, i));
      }
      p[i] = tmp[to_int(GETV(x, i))];
      
    }
  }

  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rbnbinom(
    const int& n,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta
  ) {

  NumericVector x(n);
  
  bool throw_warning = false;

  for (int i = 0; i < n; i++)
    x[i] = rng_bnbinom(GETV(size, i), GETV(alpha, i), GETV(beta, i), throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

