#include <Rcpp.h>
using namespace Rcpp;

/*
 * Inverse chi-squared distribution
 * 
 * f(x) = ((tau^2 * nu/2)^(nu/2))/gamma(nu/2) * exp(-(tau^2 * nu)/(2*x)) / (x^(1+nu/2))
 * F(x) = gamma(nu/2, (tau^2 * nu)/(2*x)) / gamma(nu/2)
 * 
 */


double pdf_scinvchisq(double x, double nu, double tau) {
  double tausq = std::pow(tau, 2);
  if (nu <= 0 || tausq <= 0)
    return NAN;
  if (x > 0)
    return std::pow(tausq*nu/2, nu/2)/R::gammafn(nu/2) *
           std::exp((-tausq*nu)/(2*x)) / std::pow(x, 1+nu/2);
  else
    return 0;
}


// [[Rcpp::export]]
NumericVector cpp_dsinvchisq(NumericVector x,
                             NumericVector nu, NumericVector tau,
                             bool log_prob = false) {
  
  int n = x.length();
  int nn = nu.length();
  int nt = tau.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn, nt));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_scinvchisq(x[i % n], nu[i % nn], tau[i % nt]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);
  
  return p;
}
  
  