#include <Rcpp.h>

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
*  Inverse-Gamma distribution
*
*  Values:
*  x
*
*  Parameters:
*  alpha > 0
*  beta > 0
*
*  f(k) = (x^(-alpha-1) * exp(-1/(beta*x))) / (Gamma(alpha)*beta^alpha)
*  F(x) = gamma(alpha, 1/(beta*x)) / Gamma(alpha)
*
*  V. Witkovsky (2001) Computing the distribution of a linear
*  combination of inverted gamma variables, Kybernetika 37(1), 79-90
*
*/

double pdf_invgamma(double x, double alpha, double beta) {
  if (ISNAN(x) || ISNAN(alpha) || ISNAN(beta))
    return NA_REAL;
  if (alpha <= 0.0 || beta <= 0.0) {
    Rcpp::warning("NaNs produced");
    return NAN;
  }
  if (x <= 0.0)
    return 0.0;
  return (pow(x, -alpha-1.0) * exp(-1.0/(beta*x))) / (R::gammafn(alpha) * pow(beta, alpha));
}


// [[Rcpp::export]]
NumericVector cpp_dinvgamma(
    const NumericVector& x,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& log_prob = false
  ) {

  std::vector<int> dims;
  dims.push_back(x.length());
  dims.push_back(alpha.length());
  dims.push_back(beta.length());
  int Nmax = *std::max_element(dims.begin(), dims.end());
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_invgamma(x[i % dims[0]], alpha[i % dims[1]], beta[i % dims[2]]);

  if (log_prob)
    p = Rcpp::log(p);

  return p;
}

