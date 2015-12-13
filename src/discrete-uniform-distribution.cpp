#include <Rcpp.h>
using namespace Rcpp;

/*
 * Discrete uniform distribution
 * 
 * Values:
 * a <= x <= b
 * 
 * f(x) = 1/(b-a+1)
 * F(x) = (floor(x)-a+1)/b-a+1
 *  
 */


double pmf_dunif(double x, int min, int max) {
  if (x < min || x > max || std::floor(x) != x)
    return 0;
  return 1/(max-min+1);
}


double cdf_dunif(double x, int min, int max) {
  if (x < min || x > max)
    return 0;
  return (std::floor(x)-min+1)/(max-min+1);
}

int invcdf_dunif(double p, int min, int max) {
  if (p <= 0 || p > 1)
    return NAN;
  return std::ceil( p*(max-min+1)+min-1 );
}

int rng_dunif(int min, int max) {
  return std::ceil( R::runif(min, max) );
}


// [[Rcpp::export]]
NumericVector cpp_ddunif(NumericVector x,
                         IntegerVector min, IntegerVector max,
                         bool log_prob = false) {
  
  int n  = x.length();
  int na = min.length();
  int nb = max.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pmf_dunif(x[i % n], min[i % na], max[i % nb]);
  
  if (!log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::exp(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pdunif(NumericVector x,
                         IntegerVector min, IntegerVector max,
                         bool lower_tail = true, bool log_prob = false) {
  
  int n  = x.length();
  int na = min.length();
  int nb = max.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_dunif(x[i % n], min[i % na], max[i % nb]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
IntegerVector cpp_qdunif(NumericVector p,
                         IntegerVector min, IntegerVector max,
                         bool lower_tail = true, bool log_prob = false) {
  
  int n  = p.length();
  int na = min.length();
  int nb = max.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, na, nb));
  IntegerVector q(Nmax);
  
  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = std::exp(p[i]);
  
  if (!lower_tail)
    for (int i = 0; i < n; i++)
      p[i] = 1-p[i];
  
  for (int i = 0; i < Nmax; i++)
    q[i] = invcdf_dunif(p[i % n], min[i % na], max[i % nb]);
  
  return q;
}


// [[Rcpp::export]]
IntegerVector cpp_rdunif(int n, IntegerVector min, IntegerVector max) {
  
  double u;
  int na = min.length();
  int nb = max.length();
  IntegerVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_dunif(min[i % na], max[i % nb]);
  
  return x;
}

