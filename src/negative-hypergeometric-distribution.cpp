#include <Rcpp.h>
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


// [[Rcpp::export]]
std::vector<double> nhyper_table(
    double n, double m, double r,
    bool cumulative = false
  ) {
  
  double j, N, start_eps;
  int ni = static_cast<int>(n);
  N = m+n;
  
  std::vector<double> t(ni), h(ni), p(ni+1);
  start_eps = 1e-100;
  h[0] = start_eps * r*n/(N-r);
  t[0] = start_eps + h[0];

  for (int i = 1; i <= ni-1; i++) {
    j = static_cast<double>(i) + r;
    h[i] = h[i-1] * j*(n+r-j)/(N-j)/(j+1.0-r);
    t[i] = t[i-1] + h[i];
  }
  
  p[0] = start_eps / t[ni-1];
  
  if (cumulative) {
    for (int i = 1; i < ni; i++)
      p[i] = t[i-1] / t[ni-1];
    p[ni] = 1.0;
  } else {
    for (int i = 1; i <= ni; i++)
      p[i] = h[i-1] / t[ni-1];
  }
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_dnhyper(
    const NumericVector& x,
    const NumericVector& n,
    const NumericVector& m,
    const NumericVector& r,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(n.length());
  dims.push_back(m.length());
  dims.push_back(r.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  std::map<std::tuple<int, int, int>, std::vector<double>> memo;
  
  for (int i = 0; i < Nmax; i++) {
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    
    if (ISNAN(x[i % dims[0]]) || ISNAN(n[i % dims[1]]) ||
        ISNAN(m[i % dims[2]]) || ISNAN(r[i % dims[3]])) {
      p[i] = NA_REAL;
    } else if (r[i % dims[3]] > m[i % dims[2]] || n[i % dims[1]] < 0.0 ||
               m[i % dims[2]] < 0.0 || r[i % dims[3]] < 0.0) {
      Rcpp::warning("NaNs produced");
      p[i] = NAN;
    } else if (x[i % dims[0]] < r[i % dims[3]] ||
               x[i % dims[0]] > (n[i % dims[1]] + r[i % dims[3]]) ||
               !isInteger(x[i % dims[0]]) || !isInteger(n[i % dims[1]]) ||
               !isInteger(m[i % dims[2]]) || !isInteger(r[i % dims[3]])) {
      p[i] = 0.0;
    } else {
      
      std::vector<double>& tmp = memo[std::make_tuple(i % dims[1], i % dims[2], i % dims[3])];
      if (!tmp.size()) {
        tmp = nhyper_table(n[i % dims[1]], m[i % dims[2]], r[i % dims[3]], false);
      }
      p[i] = tmp[static_cast<int>( x[i % dims[0]] - r[i % dims[3]] )];
      
    }
  } 
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pnhyper(
    const NumericVector& x,
    const NumericVector& n,
    const NumericVector& m,
    const NumericVector& r,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(n.length());
  dims.push_back(m.length());
  dims.push_back(r.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);
  
  std::map<std::tuple<int, int, int>, std::vector<double>> memo;
  
  for (int i = 0; i < Nmax; i++) {
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    
    if (ISNAN(x[i % dims[0]]) || ISNAN(n[i % dims[1]]) ||
        ISNAN(m[i % dims[2]]) || ISNAN(r[i % dims[3]])) {
      p[i] = NA_REAL;
    } else if (r[i % dims[3]] > m[i % dims[2]] || n[i % dims[1]] < 0.0 ||
               m[i % dims[2]] < 0.0 || r[i % dims[3]] < 0.0) {
               Rcpp::warning("NaNs produced");
      p[i] = NAN;
    } else if (x[i % dims[0]] < r[i % dims[3]]) {
      p[i] = 0.0;
    } else if (x[i % dims[0]] >= (n[i % dims[1]] + r[i % dims[3]])) {
      p[i] = 1.0;
    } else {
      
      std::vector<double>& tmp = memo[std::make_tuple(i % dims[1], i % dims[2], i % dims[3])];
      if (!tmp.size()) {
        tmp = nhyper_table(n[i % dims[1]], m[i % dims[2]], r[i % dims[3]], true);
      }
      p[i] = tmp[static_cast<int>( x[i % dims[0]] - r[i % dims[3]] )];
      
    }
  } 
  
  if (!lower_tail)
    p = 1.0 - p;
  
  if (log_prob)
    p = Rcpp::log(p);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qnhyper(
    const NumericVector& p,
    const NumericVector& n,
    const NumericVector& m,
    const NumericVector& r,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  std::vector<int> dims;
  dims.push_back(p.length());
  dims.push_back(n.length());
  dims.push_back(m.length());
  dims.push_back(r.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  std::map<std::tuple<int, int, int>, std::vector<double>> memo;
  
  for (int i = 0; i < Nmax; i++) {
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    
    if (ISNAN(p[i % dims[0]]) || ISNAN(n[i % dims[1]]) ||
        ISNAN(m[i % dims[2]]) || ISNAN(r[i % dims[3]])) {
      x[i] = NA_REAL;
    } else if (p[i % dims[0]] < 0.0 || p[i % dims[0]] > 1.0 ||
               r[i % dims[3]] > m[i % dims[2]] || n[i % dims[1]] < 0.0 ||
               m[i % dims[2]] < 0.0 || r[i % dims[3]] < 0.0) {
               Rcpp::warning("NaNs produced");
      x[i] = NAN;
    } else {
      
      std::vector<double>& tmp = memo[std::make_tuple(i % dims[1], i % dims[2], i % dims[3])];
      if (!tmp.size()) {
        tmp = nhyper_table(n[i % dims[1]], m[i % dims[2]], r[i % dims[3]], true);
      }
      
      for (int j = 0; j <= static_cast<int>( n[i % dims[1]] ); j++) {
        if (tmp[j] >= p[i % dims[0]]) {
          x[i] = static_cast<double>(j) + r[i % dims[3]];
          break;
        }
      }
      
    }
  } 
  
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rnhyper(
    const int& nn,
    const NumericVector& n,
    const NumericVector& m,
    const NumericVector& r
  ) {
  
  double u;
  std::vector<int> dims;
  dims.push_back(n.length());
  dims.push_back(m.length());
  dims.push_back(r.length());
  NumericVector x(nn);
  
  std::map<std::tuple<int, int, int>, std::vector<double>> memo;
  
  for (int i = 0; i < nn; i++) {
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    
    if (ISNAN(n[i % dims[0]]) || ISNAN(m[i % dims[1]]) || ISNAN(r[i % dims[2]])) {
      x[i] = NA_REAL;
    } else if (r[i % dims[2]] > m[i % dims[1]] || n[i % dims[0]] < 0.0 ||
               m[i % dims[1]] < 0.0 || r[i % dims[2]] < 0.0) {
      Rcpp::warning("NaNs produced");
      x[i] = NAN;
    } else {
      
      std::vector<double>& tmp = memo[std::make_tuple(i % dims[0], i % dims[1], i % dims[2])];
      if (!tmp.size()) {
        tmp = nhyper_table(n[i % dims[0]], m[i % dims[1]], r[i % dims[2]], true);
      }
      
      u = rng_unif();
      
      for (int j = 0; j <= static_cast<int>( n[i % dims[0]] ); j++) {
        if (tmp[j] >= u) {
          x[i] = static_cast<double>(j) + r[i % dims[2]];
          break;
        }
      }
      
    }
  } 
  
  return x;
}

