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
*  Gamma-Poisson distribution
*
*  Values:
*  x >= 0
*
*  Parameters:
*  alpha > 0
*  beta > 0
*
*/

inline double logpmf_gpois(double x, double alpha, double beta,
                           bool& throw_warning) {
  if (ISNAN(x) || ISNAN(alpha) || ISNAN(beta))
    return x+alpha+beta;
  if (alpha <= 0.0 || beta <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (!isInteger(x) || x < 0.0 || !R_FINITE(x))
    return R_NegInf;
  double p = beta/(1.0+beta);
  return R::lgammafn(alpha+x) - (lfactorial(x) + R::lgammafn(alpha)) +
    log(p)*x + log(1.0-p)*alpha;
}

inline std::vector<double> cdf_gpois_table(double x, double alpha, double beta) {
  
  if (x < 0.0 || !R_FINITE(x) || alpha < 0.0 || beta < 0.0)
    Rcpp::stop("inadmissible values");
  
  x = floor(x);
  std::vector<double> p_tab(TO_INT(x)+1);
  double p, qa, ga, gax, xf, px, lp;
  
  p = beta/(1.0+beta);
  qa = log(pow(1.0 - p, alpha));
  ga = R::lgammafn(alpha);
  lp = log(p);
  
  // x = 0
  
  gax = ga;
  xf = 0.0;
  px = 0.0;
  p_tab[0] = exp(qa);
  
  if (x < 1.0)
    return p_tab;
  
  // x < 2
  
  gax += log(alpha);
  px += lp;
  p_tab[1] = p_tab[0] + exp(gax - ga + px + qa);
  
  if (x < 2.0)
    return p_tab;
  
  // x >= 2
  
  for (double j = 2.0; j <= x; j += 1.0) {
    gax += log(j + alpha - 1.0);
    xf += log(j);
    px += lp;
    p_tab[TO_INT(j)] = p_tab[TO_INT(j)-1] + exp(gax - (xf + ga) + px + qa);
  }
  
  return p_tab;
}

inline double rng_gpois(double alpha, double beta,
                        bool& throw_warning) {
  if (ISNAN(alpha) || ISNAN(beta) || alpha <= 0.0 || beta <= 0.0) {
    throw_warning = true;
    return NA_REAL;
  }
  double lambda = R::rgamma(alpha, beta);
  return R::rpois(lambda);
}


// [[Rcpp::export]]
NumericVector cpp_dgpois(
    const NumericVector& x,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(alpha.length());
  dims.push_back(beta.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = logpmf_gpois(GETV(x, i), GETV(alpha, i),
                        GETV(beta, i), throw_warning);

  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pgpois(
    const NumericVector& x,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(alpha.length());
  dims.push_back(beta.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  std::map<std::tuple<int, int>, std::vector<double>> memo;
  double mx = finite_max(x);
  
  for (int i = 0; i < Nmax; i++) {
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    if (ISNAN(GETV(x, i)) || ISNAN(GETV(alpha, i)) || ISNAN(GETV(beta, i))) {
      p[i] = GETV(x, i) + GETV(alpha, i) + GETV(beta, i);
    } else if (GETV(alpha, i) <= 0.0 || GETV(beta, i) <= 0.0) {
      throw_warning = true;
      p[i] = NAN;
    } else if (GETV(x, i) < 0.0) {
      p[i] = 0.0;
    } else if (GETV(x, i) == R_PosInf) {
      p[i] = 1.0;
    } else {
      
      std::vector<double>& tmp = memo[std::make_tuple(i % dims[1], i % dims[2])];
      if (!tmp.size()) {
        tmp = cdf_gpois_table(mx, GETV(alpha, i), GETV(beta, i));
      }
      p[i] = tmp[TO_INT(GETV(x, i))];
      
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
NumericVector cpp_rgpois(
    const int& n,
    const NumericVector& alpha,
    const NumericVector& beta
  ) {

  std::vector<int> dims;
  dims.push_back(alpha.length());
  dims.push_back(beta.length());
  NumericVector x(n);
  
  bool throw_warning = false;

  for (int i = 0; i < n; i++)
    x[i] = rng_gpois(GETV(alpha, i), GETV(beta, i),
                     throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

