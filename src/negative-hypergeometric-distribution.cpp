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
    bool log_prob = false
) {
  
  int nx  = x.length();
  int nr = r.length();
  int nl = n.length();
  int nm = m.length();
  int nmax = Rcpp::max(IntegerVector::create(nx, nr, nl, nm));
  NumericVector p(nmax);
  
  for (int i = 0; i < nmax; i++)
    p[i] = pmf_nhyper(x[i % nx], m[i % nm], n[i % nl], r[i % nr]);

  if (log_prob)
    for (int i = 0; i < nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pnhyper(
    const NumericVector& x,
    const NumericVector& n,
    const NumericVector& m,
    const NumericVector& r,
    bool lower_tail = true, bool log_prob = false
) {
  
  int nx  = x.length();
  int nr = r.length();
  int nl = n.length();
  int nm = m.length();
  int nmax = Rcpp::max(IntegerVector::create(nx, nr, nl, nm));
  NumericVector p(nmax);
  
  for (int i = 0; i < nmax; i++)
    p[i] = cdf_nhyper(x[i % nx], m[i % nm], n[i % nl], r[i % nr]);
  
  if (!lower_tail)
    for (int i = 0; i < nmax; i++)
      p[i] = 1.0 - p[i];
  
  if (log_prob)
    for (int i = 0; i < nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qnhyper(
    const NumericVector& p,
    const NumericVector& n,
    const NumericVector& m,
    const NumericVector& r,
    bool lower_tail = true, bool log_prob = false
) {
  
  int np  = p.length();
  int nr = r.length();
  int nl = n.length();
  int nm = m.length();
  int nmax = Rcpp::max(IntegerVector::create(np, nr, nl, nm));
  NumericVector q(nmax);
  NumericVector pp = Rcpp::clone(p);
  
  if (log_prob)
    for (int i = 0; i < np; i++)
      pp[i] = exp(pp[i]);
  
  if (!lower_tail)
    for (int i = 0; i < np; i++)
      pp[i] = 1.0 - pp[i];
  
  for (int i = 0; i < nmax; i++)
    q[i] = invcdf_nhyper(pp[i % np], m[i % nm], n[i % nl], r[i % nr]);
  
  return q;
}


// [[Rcpp::export]]
NumericVector cpp_rnhyper(
    const int nn,
    const NumericVector& n,
    const NumericVector& m,
    const NumericVector& r
) {
  
  double u;
  int nr = r.length();
  int nl = n.length();
  int nm = m.length();
  NumericVector x(nn);
  
  /*
  if (nr == 1 && nl == 1 && nm == 1) {
    
    if (ISNAN(r[0]) || ISNAN(n[0]) || ISNAN(m[0])) {
      for (int j = 0; j < nn; j++)
        x[j] = NA_REAL;
      return x;
    }
    
    double p1, rmax, h, t, i;
    double N = m[0]+n[0];
    
    h = 1.0e-100;
    rmax = n[0]+r[0];
    t = h;
    
    if (r[0] > rmax || m[0] > N) {
      Rcpp::warning("NanN produced");
      return NAN;
    }
    
    if (r[0] > rmax) {
      for (int j = 0; j < nn; j++)
        x[j] = NAN;
      return x;
    }

    i = r[0];
    do {
      h *= i*(n[0]+r[0]-i)/(N-i)/(i+1.0-r[0]);
      t += h;
      i += 1.0;
    } while (i <= rmax-1.0);
    
    h = 1.0e-100;
    p1 = h;
    
    std::vector<double> ptab(static_cast<int>(rmax)+1, 0.0);
    
    i = r[0];
    do {
      h *= i*(n[0]+r[0]-i)/(N-i)/(i+1.0-r[0]);
      p1 += h;
      i += 1.0;
      ptab[static_cast<int>(i)] = p1/t;
    } while (i <= rmax);
    
    for (int j = 0; j < nn; j++) {
      u = rng_unif();
      
      for (int z = rmax; z >= 0; z--) {
        if (u > ptab[z]) {
          x[j] = static_cast<double>(z+1);
          break;
        }
      }
      
    }
    
  } else {
    */
    
    for (int i = 0; i < nn; i++) {
      u = rng_unif();
      x[i] = invcdf_nhyper(u, m[i % nm], n[i % nl], r[i % nr]);
    }
    
  //}
  
  return x;
}

