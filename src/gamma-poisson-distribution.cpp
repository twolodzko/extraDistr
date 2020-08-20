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
#ifdef IEEE_754
  if (ISNAN(x) || ISNAN(alpha) || ISNAN(beta))
    return x+alpha+beta;
#endif
  if (alpha <= 0.0 || beta <= 0.0) {
    throw_warning = true;
    return NAN;
  }
  if (!isInteger(x) || x < 0.0 || !R_FINITE(x))
    return R_NegInf;
  // p = beta/(1.0+beta);
  double p = exp( log(beta) - log1p(beta) );
  return R::lgammafn(alpha+x) - lfactorial(x) - R::lgammafn(alpha) +
    log(p)*x + log(1.0-p)*alpha;
}

inline std::vector<double> cdf_gpois_table(double x, double alpha, double beta) {
  
  if (x < 0.0 || !R_FINITE(x) || alpha < 0.0 || beta < 0.0)
    Rcpp::stop("inadmissible values");
  
  int ix = to_pos_int(x);
  std::vector<double> p_tab(ix+1);
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
  
  if (ix < 1)
    return p_tab;
  
  // x < 2
  
  gax += log(alpha);
  px += lp;
  p_tab[1] = p_tab[0] + exp(gax - ga + px + qa);
  
  if (ix < 2)
    return p_tab;
  
  // x >= 2
  
  double dj;
  
  for (int j = 2; j <= ix; j++) {
    if (j % 10000 == 0)
      Rcpp::checkUserInterrupt();
    dj = to_dbl(j);
    gax += log(dj + alpha - 1.0);
    xf += log(dj);
    px += lp;
    p_tab[j] = p_tab[j-1] + exp(gax - (xf + ga) + px + qa);
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
  
  if (std::min({x.length(), alpha.length(), beta.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    alpha.length(),
    beta.length()
  });
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
  
  if (std::min({x.length(), alpha.length(), beta.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    alpha.length(),
    beta.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  std::map<std::tuple<int, int>, std::vector<double>> memo;
  double mx = finite_max_int(x);
  
  for (int i = 0; i < Nmax; i++) {
    if (i % 100 == 0)
      Rcpp::checkUserInterrupt();
    
#ifdef IEEE_754
    if (ISNAN(GETV(x, i)) || ISNAN(GETV(alpha, i)) || ISNAN(GETV(beta, i))) {
      p[i] = GETV(x, i) + GETV(alpha, i) + GETV(beta, i);
      continue;
    }
#endif
    
    if (GETV(alpha, i) <= 0.0 || GETV(beta, i) <= 0.0) {
      throw_warning = true;
      p[i] = NAN;
    } else if (GETV(x, i) < 0.0) {
      p[i] = 0.0;
    } else if (GETV(x, i) == R_PosInf) {
      p[i] = 1.0;
    } else if (is_large_int(GETV(x, i))) {
      p[i] = NA_REAL;
      Rcpp::warning("NAs introduced by coercion to integer range");
    } else {
      
      std::vector<double>& tmp = memo[std::make_tuple(
        static_cast<int>(i % alpha.length()),
        static_cast<int>(i % beta.length())
      )];
      
      if (!tmp.size()) {
        tmp = cdf_gpois_table(mx, GETV(alpha, i), GETV(beta, i));
      }
      p[i] = tmp[to_pos_int(GETV(x, i))];
      
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
  
  if (std::min({alpha.length(), beta.length()}) < 1) {
    Rcpp::warning("NAs produced");
    return NumericVector(n, NA_REAL);
  }

  NumericVector x(n);
  
  bool throw_warning = false;

  for (int i = 0; i < n; i++)
    x[i] = rng_gpois(GETV(alpha, i), GETV(beta, i),
                     throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

