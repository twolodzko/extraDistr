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
using std::sin;
using std::cos;
using std::tan;
using std::atan;
using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;


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

double pmf_bnbinom(double k, double r, double alpha, double beta) {
  if (ISNAN(k) || ISNAN(r) || ISNAN(alpha) || ISNAN(beta))
    return NA_REAL;
  if (alpha <= 0.0 || beta <= 0.0 || r < 0.0 || !isInteger(r, false)) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(k) || k < 0.0 || !R_FINITE(k))
    return 0.0;
  return (R::gammafn(r+k) / (R::gammafn(k+1.0) * R::gammafn(r))) *
          R::beta(alpha+r, beta+k) / R::beta(alpha, beta);
}

double logpmf_bnbinom(double k, double r, double alpha, double beta) {
  if (ISNAN(k) || ISNAN(r) || ISNAN(alpha) || ISNAN(beta))
    return NA_REAL;
  if (alpha <= 0.0 || beta <= 0.0 || r < 0.0 || !isInteger(r, false)) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (!isInteger(k) || k < 0.0 || !R_FINITE(k))
    return R_NegInf;
  return (R::lgammafn(r+k) - (R::lgammafn(k+1.0) + R::lgammafn(r))) +
    R::lbeta(alpha+r, beta+k) - R::lbeta(alpha, beta);
}

std::vector<double> cdf_bnbinom_table(double k, double r, double alpha, double beta) {

  k = floor(k);
  std::vector<double> p_tab(static_cast<int>(k)+1);
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
  
  if (k < 1.0)
    return p_tab;
  
  // k < 2
  
  grx += log(r);
  gbx += log(beta);
  gabrx += log(alpha + beta + r);
  p_tab[1] = p_tab[0] + exp(grx - gr + gar + gbx - gabrx - bab);
  
  if (k < 2.0)
    return p_tab;
  
  // k >= 2
  
  for (double j = 2.0; j <= k; j += 1.0) {
    grx += log(r + j - 1.0);
    gbx += log(beta + j - 1.0);
    gabrx += log(alpha + beta + r + j - 1.0);
    xf += log(j);
    p_tab[static_cast<int>(j)] = p_tab[static_cast<int>(j)-1] +
      exp(grx - (xf + gr) + gar + gbx - gabrx - bab);
  }
  
  return p_tab;
}

double rng_bnbinom(double r, double alpha, double beta) {
  if (ISNAN(r) || ISNAN(alpha) || ISNAN(beta) ||
      alpha <= 0.0 || beta <= 0.0 || r < 0.0 || !isInteger(r, false)) {
    Rcpp::warning("NAs produced");
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

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(size.length());
  dims.push_back(alpha.length());
  dims.push_back(beta.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = logpmf_bnbinom(x[i % dims[0]], size[i % dims[1]],
                          alpha[i % dims[2]], beta[i % dims[3]]);

  if (!log_prob)
    p = Rcpp::exp(p);

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

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(size.length());
  dims.push_back(alpha.length());
  dims.push_back(beta.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);

  std::map<std::tuple<int, int, int>, std::vector<double>> memo;
  double mx = finite_max(x);
  
  for (int i = 0; i < Nmax; i++) {
    
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    
    if (ISNAN(x[i % dims[0]]) || ISNAN(size[i % dims[1]]) ||
        ISNAN(alpha[i % dims[2]]) || ISNAN(beta[i % dims[3]])) {
      p[i] = NA_REAL;
    } else if (alpha[i % dims[2]] <= 0.0 || beta[i % dims[3]] <= 0.0 ||
               size[i % dims[1]] < 0.0 || !isInteger(size[i % dims[1]], false)) {
      Rcpp::warning("NaNs produced");
      p[i] = NAN;
    } else if (x[i % dims[0]] < 0.0) {
      p[i] = 0.0;
    } else if (!R_FINITE(x[i % dims[0]])) {
      p[i] = 1.0;
    } else {
      
      std::vector<double>& tmp = memo[std::make_tuple(i % dims[1], i % dims[2], i % dims[3])];
      if (!tmp.size()) {
        tmp = cdf_bnbinom_table(mx, size[i % dims[1]],
                                alpha[i % dims[2]], beta[i % dims[3]]);
      }
      p[i] = tmp[static_cast<int>(x[i % dims[0]])];
      
    }
  }

  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rbnbinom(
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

  for (int i = 0; i < n; i++)
    x[i] = rng_bnbinom(size[i % dims[0]], alpha[i % dims[1]], beta[i % dims[2]]);

  return x;
}

