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



double pmf_nhyper(double x, double m, double n, double r) {
  
  if (ISNAN(x) || ISNAN(r) || ISNAN(n) || ISNAN(m))
    return NA_REAL;
  
  double p, rmax, h, t, i;
  double N = m+n;
  
  h = 1e-100;
  rmax = n+r;
  p = h;
  t = h;
  
  if (r > rmax || m > N) {
    Rcpp::warning("NanN produced");
    return NAN;
  }
  
  if (x < r || x > rmax ||
      !isInteger(x) || !isInteger(r) ||
      !isInteger(n) || !isInteger(m))
      return 0.0;
  
  if (!tol_equal(x, r)) {
    i = r;
    do {
      h *= i*(n+r-i)/(N-i)/(i+1.0-r);
      t += h;
      i += 1.0;
    } while (i <= x-1.0);
    p = h;
    
    if (tol_equal(x, rmax))
      return p/t;
  }
  
  i = x;
  do {
    h *= i*(n+r-i)/(N-i)/(i+1.0-r);
    t += h;
    i += 1.0;
  } while (i <= rmax-1.0);
  
  return p/t;
}


double cdf_nhyper(double x, double m, double n, double r) {
  
  if (ISNAN(x) || ISNAN(r) || ISNAN(n) || ISNAN(m))
    return NA_REAL;
  
  double p, rmax, h, t, i;
  double N = m+n;
  
  h = 1e-100;
  rmax = n+r;
  p = h;
  
  if (r > rmax || m > N) {
    Rcpp::warning("NanN produced");
    return NAN;
  }
  
  if (x < r)
    return 0.0;
  
  if (x >= rmax)
    return 1.0;
  
  if (tol_equal(x, r)) {
    t = h;
  } else {
    i = r;
    do {
      h *= i*(n+r-i)/(N-i)/(i+1.0-r);
      p += h;
      i += 1.0;
    } while (i <= x-1.0);
    t = p;
    
    if (tol_equal(x, rmax))
      return p/t;
  }
  
  i = x;
  do {
    h *= i*(n+r-i)/(N-i)/(i+1.0-r);
    t += h;
    i += 1.0;
  } while (i <= rmax-1.0);
  
  return p/t;
}


double invcdf_nhyper(double p, double m, double n, double r) {
  
  if (ISNAN(p) || ISNAN(r) || ISNAN(n) || ISNAN(m))
    return NA_REAL;
  
  double p1, rmax, h, t, i, pt;
  double N = m+n;
  rmax = n+r;
  
  if (r > rmax || m > N) {
    Rcpp::warning("NanN produced");
    return NAN;
  }
  
  h = 1e-100;
  t = h;
  
  i = r;
  do {
    h *= i*(n+r-i)/(N-i)/(i+1.0-r);
    t += h;
    i += 1.0;
  } while (i <= rmax-1.0);
  
  pt = p*t;
  h = 1e-100;
  p1 = h;
  
  // Rcpp::Rcout << "** I = " << i << " ; P = " << p1/t << std::endl;
  if (p1 >= pt)
    return r;
  
  i = r;
  do {
    h *= i*(n+r-i)/(N-i)/(i+1.0-r);
    p1 += h;
    i += 1.0;
    // Rcpp::Rcout << "I = " << i << " ; P = " << p1/t << std::endl;
    if (p1 >= pt)
      break;
  } while (i <= rmax);
  
  return i;
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
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pmf_nhyper(x[i % dims[0]], m[i % dims[1]], n[i % dims[2]], r[i % dims[3]]);

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
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_nhyper(x[i % dims[0]], m[i % dims[1]], n[i % dims[2]], r[i % dims[3]]);
  
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
  
  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_nhyper(pp[i % dims[0]], m[i % dims[1]], n[i % dims[2]], r[i % dims[3]]);
  
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
  
  for (int i = 0; i < nn; i++) {
    u = rng_unif();
    x[i] = invcdf_nhyper(u, m[i % dims[0]], n[i % dims[1]], r[i % dims[2]]);
  }

  return x;
}

