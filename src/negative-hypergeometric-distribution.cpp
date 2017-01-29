#include <Rcpp.h>
#include "shared.h"
// [[Rcpp::plugins(cpp11)]]

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using Rcpp::NumericVector;


std::vector<double> nhyper_table(
    double n, double m, double r,
    bool cumulative = false
  ) {
  
  if (n < 0.0 || m < 0.0 || r < 0.0 || r > m)
    Rcpp::stop("inadmissible values");
  
  double j, N, start_eps;
  long int ni = to_int(n);
  N = m+n;
  
  std::vector<double> t(ni), h(ni), p(ni+1);
  start_eps = 1e-100;
  h[0] = start_eps * r*n/(N-r);
  t[0] = start_eps + h[0];

  for (long int i = 1; i <= ni-1; i++) {
    j = to_dbl(i) + r;
    h[i] = h[i-1] * j*(n+r-j)/(N-j)/(j+1.0-r);
    t[i] = t[i-1] + h[i];
  }
  
  p[0] = start_eps / t[ni-1];
  
  if (cumulative) {
    for (long int i = 1; i < ni; i++)
      p[i] = t[i-1] / t[ni-1];
    p[ni] = 1.0;
  } else {
    for (long int i = 1; i <= ni; i++)
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
  
  int Nmax = std::max({
    x.length(),
    n.length(),
    m.length(),
    r.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  std::map<std::tuple<int, int, int>, std::vector<double>> memo;
  
  for (int i = 0; i < Nmax; i++) {
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    
    if (ISNAN(GETV(x, i)) || ISNAN(GETV(n, i)) ||
        ISNAN(GETV(m, i)) || ISNAN(GETV(r, i))) {
      p[i] = GETV(x, i) + GETV(n, i) + GETV(m, i) + GETV(r, i);
    } else if (GETV(r, i) > GETV(m, i) || GETV(n, i) < 0.0 ||
               GETV(m, i) < 0.0 || GETV(r, i) < 0.0 ||
               !isInteger(GETV(n, i), false) ||
               !isInteger(GETV(m, i), false) ||
               !isInteger(GETV(r, i), false)) {
      throw_warning = true;
      p[i] = NAN;
    } else if (!isInteger(GETV(x, i)) || GETV(x, i) < GETV(r, i) ||
               GETV(x, i) > (GETV(n, i) + GETV(r, i))) {
      p[i] = 0.0;
    } else {
      
      std::vector<double>& tmp = memo[std::make_tuple(i % n.length(),
                                                      i % m.length(),
                                                      i % r.length())];
      if (!tmp.size()) {
        tmp = nhyper_table(GETV(n, i), GETV(m, i), GETV(r, i), false);
      }
      p[i] = tmp[to_int( GETV(x, i) - GETV(r, i) )];
      
    }
  } 
  
  if (log_prob)
    p = Rcpp::log(p);
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
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
  
  int Nmax = std::max({
    x.length(),
    n.length(),
    m.length(),
    r.length()
  });
  NumericVector p(Nmax);
  
  bool throw_warning = false;
  
  std::map<std::tuple<int, int, int>, std::vector<double>> memo;
  
  for (int i = 0; i < Nmax; i++) {
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    
    if (ISNAN(GETV(x, i)) || ISNAN(GETV(n, i)) ||
        ISNAN(GETV(m, i)) || ISNAN(GETV(r, i))) {
      p[i] = GETV(x, i) + GETV(n, i) + GETV(m, i) + GETV(r, i);
    } else if (GETV(r, i) > GETV(m, i) || GETV(n, i) < 0.0 ||
               GETV(m, i) < 0.0 || GETV(r, i) < 0.0 ||
               !isInteger(GETV(n, i), false) ||
               !isInteger(GETV(m, i), false) ||
               !isInteger(GETV(r, i), false)) {
      throw_warning = true;
      p[i] = NAN;
    } else if (GETV(x, i) < GETV(r, i)) {
      p[i] = 0.0;
    } else if (GETV(x, i) >= (GETV(n, i) + GETV(r, i))) {
      p[i] = 1.0;
    } else {
      
      std::vector<double>& tmp = memo[std::make_tuple(i % n.length(),
                                                      i % m.length(),
                                                      i % r.length())];
      if (!tmp.size()) {
        tmp = nhyper_table(GETV(n, i), GETV(m, i), GETV(r, i), true);
      }
      p[i] = tmp[to_int( GETV(x, i) - GETV(r, i) )];
      
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
NumericVector cpp_qnhyper(
    const NumericVector& p,
    const NumericVector& n,
    const NumericVector& m,
    const NumericVector& r,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {
  
  int Nmax = std::max({
    p.length(),
    n.length(),
    m.length(),
    r.length()
  });
  NumericVector x(Nmax);
  NumericVector pp = Rcpp::clone(p);
  
  bool throw_warning = false;
  
  if (log_prob)
    pp = Rcpp::exp(pp);
  
  if (!lower_tail)
    pp = 1.0 - pp;
  
  std::map<std::tuple<int, int, int>, std::vector<double>> memo;
  
  for (int i = 0; i < Nmax; i++) {
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    
    if (ISNAN(GETV(pp, i)) || ISNAN(GETV(n, i)) ||
        ISNAN(GETV(m, i)) || ISNAN(GETV(r, i))) {
      x[i] = GETV(pp, i) + GETV(n, i) + GETV(m, i) + GETV(r, i);
    } else if (!VALID_PROB(GETV(pp, i)) ||
               GETV(r, i) > GETV(m, i) || GETV(n, i) < 0.0 ||
               GETV(m, i) < 0.0 || GETV(r, i) < 0.0 ||
               !isInteger(GETV(n, i), false) ||
               !isInteger(GETV(m, i), false) ||
               !isInteger(GETV(r, i), false)) {
      throw_warning = true;
      x[i] = NAN;
    } else {
      
      std::vector<double>& tmp = memo[std::make_tuple(i % n.length(),
                                                      i % m.length(),
                                                      i % r.length())];
      if (!tmp.size()) {
        tmp = nhyper_table(GETV(n, i), GETV(m, i), GETV(r, i), true);
      }
      
      for (int j = 0; j <= to_int( GETV(n, i) ); j++) {
        if (tmp[j] >= GETV(pp, i)) {
          x[i] = to_dbl(j) + GETV(r, i);
          break;
        }
      }
      
    }
  } 
  
  if (throw_warning)
    Rcpp::warning("NaNs produced");
  
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
  NumericVector x(nn);
  
  bool throw_warning = false;
  
  std::map<std::tuple<int, int, int>, std::vector<double>> memo;
  
  for (int i = 0; i < nn; i++) {
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();
    
    if (ISNAN(GETV(n, i)) || ISNAN(GETV(m, i)) || ISNAN(GETV(r, i)) ||
        GETV(r, i) > GETV(m, i) || GETV(n, i) < 0.0 ||
        GETV(m, i) < 0.0 || GETV(r, i) < 0.0 || !isInteger(GETV(n, i), false) ||
        !isInteger(GETV(m, i), false) || !isInteger(GETV(r, i), false)) {
      throw_warning = true;
      x[i] = NA_REAL;
    } else {
      
      std::vector<double>& tmp = memo[std::make_tuple(i % n.length(),
                                                      i % m.length(),
                                                      i % r.length())];
      if (!tmp.size()) {
        tmp = nhyper_table(GETV(n, i), GETV(m, i), GETV(r, i), true);
      }
      
      u = rng_unif();
      
      for (int j = 0; j <= to_int( GETV(n, i) ); j++) {
        if (tmp[j] >= u) {
          x[i] = to_dbl(j) + GETV(r, i);
          break;
        }
      }
      
    }
  } 
  
  if (throw_warning)
    Rcpp::warning("NAs produced");
  
  return x;
}

