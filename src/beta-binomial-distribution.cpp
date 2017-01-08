#include <Rcpp.h>
#include "const.h"
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
*  Beta-binomial distribution
*
*  Values:
*  x
*
*  Parameters:
*  k > 0
*  alpha > 0
*  beta > 0
*
*  f(k) = choose(n, k) * (beta(k+alpha, n-k+beta)) / (beta(alpha, beta))
*
*/

double pmf_bbinom(double k, double n, double alpha,
                  double beta, bool& throw_warning) {
  if (ISNAN(k) || ISNAN(n) || ISNAN(alpha) || ISNAN(beta))
    return k+n+alpha+beta;
  if (alpha < 0.0 || beta < 0.0 || n < 0.0 || !isInteger(n, false)) {
    throw_warning = true;
    return NAN;
  }
  if (!isInteger(k) || k < 0.0 || k > n)
    return 0.0;
  return R::choose(n, k) * R::beta(k+alpha, n-k+beta) / R::beta(alpha, beta);
}

double logpmf_bbinom(double k, double n, double alpha,
                     double beta, bool& throw_warning) {
  if (ISNAN(k) || ISNAN(n) || ISNAN(alpha) || ISNAN(beta))
    return k+n+alpha+beta;
  if (alpha < 0.0 || beta < 0.0 || n < 0.0 || !isInteger(n, false)) {
    throw_warning = true;
    return NAN;
  }
  if (!isInteger(k) || k < 0.0 || k > n)
    return R_NegInf;
  return R::lchoose(n, k) + R::lbeta(k+alpha, n-k+beta) - R::lbeta(alpha, beta);
}

std::vector<double> cdf_bbinom_table(double k, double n, double alpha, double beta) {
  
  if (k < 0.0 || k > n || alpha < 0.0 || beta < 0.0)
    Rcpp::stop("inadmissible values");

  k = floor(k);
  std::vector<double> p_tab(static_cast<int>(k)+1);
  double nck, bab, gx, gy, gxy;
  
  bab = R::lbeta(alpha, beta);
  gxy = R::lgammafn(alpha + beta + n);
  
  // k = 0
  
  nck = 0.0;
  gx = R::lgammafn(alpha);
  gy = R::lgammafn(beta + n);
  p_tab[0] = exp(nck + gx + gy - gxy - bab);
  
  if (k < 1.0)
    return p_tab;
  
  // k < 2
  
  nck += log(n);
  gx += log(alpha);
  gy -= log(n + beta - 1.0);
  p_tab[1] = p_tab[0] + exp(nck + gx + gy - gxy - bab);
  
  if (k < 2.0)
    return p_tab;
  
  // k >= 1
  
  for (double j = 2.0; j <= k; j += 1.0) {
    nck += log((n + 1.0 - j)/j);
    gx += log(j + alpha - 1.0);
    gy -= log(n + beta - j);
    p_tab[static_cast<int>(j)] = p_tab[static_cast<int>(j)-1] + exp(nck + gx + gy - gxy - bab);
  }
  
  return p_tab;
}

double rng_bbinom(double n, double alpha,
                  double beta, bool& throw_warning) {
  if (ISNAN(n) || ISNAN(alpha) || ISNAN(beta) ||
      alpha < 0.0 || beta < 0.0 || n < 0.0 || !isInteger(n, false)) {
    throw_warning = true;
    return NA_REAL;
  }
  double prob = R::rbeta(alpha, beta);
  return R::rbinom(n, prob);
}


// [[Rcpp::export]]
NumericVector cpp_dbbinom(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(size.length());
  dims.push_back(alpha.length());
  dims.push_back(beta.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = logpmf_bbinom(x[i % dims[0]], size[i % dims[1]],
                         alpha[i % dims[2]], beta[i % dims[3]],
                         throw_warning);

  if (!log_prob)
    p = Rcpp::exp(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
NumericVector cpp_pbbinom(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(size.length());
  dims.push_back(alpha.length());
  dims.push_back(beta.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  std::map<std::tuple<int, int, int>, std::vector<double>> memo;
  double mx = std::min(finite_max(x), finite_max(size));
  
  for (int i = 0; i < Nmax; i++) {
    
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    
    if (ISNAN(x[i % dims[0]]) || ISNAN(size[i % dims[1]]) ||
        ISNAN(alpha[i % dims[2]]) || ISNAN(beta[i % dims[3]])) {
      p[i] = x[i % dims[0]] + size[i % dims[1]] + alpha[i % dims[2]] + beta[i % dims[3]];
    } else if (alpha[i % dims[2]] <= 0.0 || beta[i % dims[3]] <= 0.0 ||
               size[i % dims[1]] < 0.0 || !isInteger(size[i % dims[1]], false)) {
      throw_warning = true;
      p[i] = NAN;
    } else if (x[i % dims[0]] < 0.0) {
      p[i] = 0.0;
    } else if (x[i % dims[0]] >= size[i % dims[1]]) {
      p[i] = 1.0;
    } else {
      
      std::vector<double>& tmp = memo[std::make_tuple(i % dims[1], i % dims[2], i % dims[3])];
      if (!tmp.size()) {
        tmp = cdf_bbinom_table(mx, size[i % dims[1]],
                               alpha[i % dims[2]], beta[i % dims[3]]);
      }
      p[i] = tmp[static_cast<int>(x[i % dims[0]])];
      
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
NumericVector cpp_rbbinom(
    const int& n,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta
  ) {

  std::vector<int> dims;
  dims.push_back(size.length());
  dims.push_back(alpha.length());
  dims.push_back(beta.length());
  NumericVector x(n);
  
  bool throw_warning = false;

  for (int i = 0; i < n; i++)
    x[i] = rng_bbinom(size[i % dims[0]], alpha[i % dims[1]], beta[i % dims[2]], throw_warning);
  
  if (throw_warning)
    Rcpp::warning("NAs produced");

  return x;
}

