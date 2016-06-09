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
* Multinomial distribution
* 
* x[i]       number of values of i-th category drawn
* n = sum(x) total number of draws
* p[i]       probability of drawing i-th category value
*  
* f(x) = n!/prod(x[i]!) * prod(p[i]^x[i])
*
*/


// [[Rcpp::export]]
NumericVector cpp_dmnom(
    const NumericMatrix& x,
    const NumericVector& size,
    const NumericMatrix& prob,
    bool log_prob = false
  ) {
  
  int n = x.nrow();
  int m = x.ncol();
  int k = prob.ncol();
  int np = prob.nrow();
  int ns = size.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, ns, np));
  NumericVector p(Nmax);
  
  if (m != k)
    Rcpp::stop("Number of columns in 'x' does not equal number of columns in 'prob'.");

  for (int i = 0; i < Nmax; i++) {
    
    double n_fac = lfactorial(size[i % ns]);
    double prod_xfac = 0;
    double prod_pow_px = 0;
    
    double sum_x = 0;
    double sum_p = 0;
    bool wrong_p = false;
    bool wrong_x = false;

    for (int j = 0; j < k; j++) {
      if (prob(i % np, j) < 0 || prob(i % np, j) > 1) {
        wrong_p = true;
        break;
      }
      if (x(i % n, j) < 0 || !isInteger(x(i % n, j))) {
        wrong_x = true;
      } else {
        sum_x += x(i % n, j);
        prod_xfac += lfactorial(x(i % n, j));
        prod_pow_px += log(prob(i % np, j)) * x(i % n, j);
      }
      sum_p += prob(i % np, j)*P_NORM_CONST;
    }

    if (!tol_equal(sum_p/P_NORM_CONST, 1.0) || wrong_p) {
      Rcpp::warning("NaNs produced");
      p[i] = NAN; 
    } else if (floor(size[i % ns]) != size[i % ns] ||
        sum_x < 0 || sum_x != size[i % ns] || wrong_x) {
      p[i] = -INFINITY;
    } else {
      p[i] = n_fac - prod_xfac + prod_pow_px;
    }
  }

  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = exp(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericMatrix cpp_rmnom(
    const int n,
    const NumericVector& size,
    const NumericMatrix& prob
  ) {
  
  int k = prob.ncol();
  int np = prob.nrow();
  int ns = size.length();
  
  NumericMatrix x(n, k);
  
  for (int i = 0; i < n; i++) {
    
    int size_left = size[i % ns];
    double sum_p = 0;
    bool wrong_p = false;
    
    for (int j = 0; j < k-1; j++) {
      if (prob(i % np, j) < 0 || prob(i % np, j) > 1) {
        wrong_p = true;
        break;
      }
      sum_p += prob(i % np, j)*P_NORM_CONST;
      x(i, j) = R::rbinom(size_left, prob(i % np, j));
      size_left -= x(i, j);
    }
    
    x(i, k-1) = size_left;
    sum_p += prob(i % np, k-1)*P_NORM_CONST;
    
    if (!tol_equal(sum_p/P_NORM_CONST, 1.0) || wrong_p) {
      Rcpp::warning("NaNs produced");
      for (int j = 0; j < k; j++)
        x(i, j) = NAN;
    }
  }
  
  return x;
}

