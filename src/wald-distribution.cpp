#include <Rcpp.h>
#include "shared.h"
using namespace Rcpp;

/*
 * Wald distribution
 * 
 * Parameters:
 * mu > 0
 * lambda > 0
 * 
 * Values:
 * x > 0
 *
 * 
 */

double pdf_wald(double x, double mu, double lambda) {
  if (mu <= 0 || lambda <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return sqrt(lambda/(2*PI*pow(x, 3))) *
         exp((-lambda*pow(x-mu, 2))/(2*pow(mu, 2)*x));
}

double cdf_wald(double x, double mu, double lambda) {
  if (mu <= 0 || lambda <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  return Phi(sqrt(lambda/x)*(x/mu-1)) +
         exp((2*lambda)/mu) *
         Phi(-sqrt(lambda/x)*(x/mu+1));
}

double rng_wald(double mu, double lambda) {
  if (mu <= 0 || lambda <= 0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  double u = R::runif(0, 1);
  double y = pow(R::rnorm(0, 1), 2);
  double x = mu + (pow(mu, 2)*y)/(2*lambda) - mu/(2*lambda) *
             sqrt(4*mu*lambda*y+pow(mu, 2)*pow(y, 2));
  if (u <= mu/(mu+x))
    return x;
  else
    return pow(mu, 2)/x;
}


// [[Rcpp::export]]
NumericVector cpp_dwald(NumericVector x,
                        NumericVector mu, NumericVector lambda,
                        bool log_prob = false) {
  
  int n  = x.length();
  int nm = mu.length();
  int nl = lambda.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, nl));
  NumericVector p(Nmax);
  double z;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_wald(x[i % n], mu[i % nm], lambda[i % nl]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pwald(NumericVector x,
                        NumericVector mu, NumericVector lambda,
                        bool lower_tail = true, bool log_prob = false) {
  
  int n  = x.length();
  int nm = mu.length();
  int nl = lambda.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nm, nl));
  NumericVector p(Nmax);
  double z;
  
  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_wald(x[i % n], mu[i % nm], lambda[i % nl]);
  
  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);
  
  return p;
}


// [[Rcpp::export]]
NumericVector cpp_rwald(int n,
                        NumericVector mu, NumericVector lambda) {
  
  int nm = mu.length();
  int nl = lambda.length();
  NumericVector x(n);
  
  for (int i = 0; i < n; i++)
    x[i] = rng_wald(mu[i % nm], lambda[i % nl]);
  
  return x;
}

