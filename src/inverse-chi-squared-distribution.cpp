#include <Rcpp.h>
using namespace Rcpp;

/*
 * Inverse chi-squared distribution
 * 
 * f(x) = (2^(-nu/2))/gamma(nu/2) * x^(-nu/2-1) * exp(-1/(2*x))
 * F(x) = gamma(nu/2, 1/(2*x)) / gamma(nu/2)
 * 
 */


double pdf_invchisq(double x, double nu) {
  if (nu <= 0)
    return NAN;
  if (x > 0)
    return std::pow(2, -nu/2) / R::gammafn(nu/2) * std::pow(x, -nu/2-1) * std::exp(-1/(2*x));
  else
    return 0;
}


// double pdf_invchisq(double x, double nu, double tau) {
//   if (nu <= 0 || tau <= 0)
//     return NAN;
//   if (x > 0)
//     return std::pow((tau*nu)/2, nu/2) * (std::exp((-nu*tau)/(2*x)) / std::pow(x, 1+nu/2));
//   else
//     return 0;
// }


// [[Rcpp::export]]
NumericVector cpp_dinvchisq(NumericVector x,
                            NumericVector nu,
                            bool log_prob = false) {
  
  int n = x.length();
  int nn = nu.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nn));
  NumericVector p(Nmax);
  
  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_invchisq(x[i % n], nu[i % nn]);
  
  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = std::log(p[i]);
  
  return p;
}
  
  