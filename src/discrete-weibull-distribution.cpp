#include <Rcpp.h>
using namespace Rcpp;


/*
*  Discrete Weibull distribution
*
*  Values:
*  x >= 0
*
*  Parameters:
*  0 < q < 1
*  beta
*
*  f(x)    = q^x^beta - q^(x+1)^beta
*  F(x)    = 1-q^(x+1)^beta
*  F^-1(p) = ceiling(pow(log(1-p)/log(q), 1/beta) - 1)
*
*  Nakagawa and Osaki (1975), "The Discrete Weibull Distribution",
*  IEEE Transactions on Reliability, R-24, pp. 300-301.
*
*/

double pdf_dweibull(int x, double q, double beta) {
  if (x >= 0) {
    return pow(q, pow(x, beta)) - pow(q, pow(x+1, beta));
  } else {
    return 0;
  }
}

double cdf_dweibull(int x, double q, double beta) {
  if (x >= 0) {
    return 1-pow(q, pow(x+1, beta));
  } else {
    return 0;
  }
}

double invcdf_dweibull(double p, double q, double beta) {
  return std::ceil(pow(log(1-p)/log(q), 1/beta) - 1);
}


// [[Rcpp::export]]
NumericVector cpp_ddweibull(IntegerVector x,
                            NumericVector q, NumericVector beta,
                            bool log_prob = false) {

  if (is_true(any(q <= 0)) || is_true(any(q >= 1)))
    throw Rcpp::exception("Values of q should be 0 < q < 1.");
  if (is_true(any(beta <= 0)))
    throw Rcpp::exception("Values of beta should be > 0.");

  int n  = x.length();
  int nq = q.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nq, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = pdf_dweibull(x[i % n], q[i % nq], beta[i % nb]);

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_pdweibull(IntegerVector x,
                            NumericVector q, NumericVector beta,
                            bool lower_tail = true, bool log_prob = false) {

  if (is_true(any(q <= 0)) || is_true(any(q >= 1)))
    throw Rcpp::exception("Values of q should be 0 < q < 1.");
  if (is_true(any(beta <= 0)))
    throw Rcpp::exception("Values of beta should be > 0.");

  int n  = x.length();
  int nq = q.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nq, nb));
  NumericVector p(Nmax);

  for (int i = 0; i < Nmax; i++)
    p[i] = cdf_dweibull(x[i % n], q[i % nq], beta[i % nb]);

  if (!lower_tail)
    for (int i = 0; i < Nmax; i++)
      p[i] = 1-p[i];

  if (log_prob)
    for (int i = 0; i < Nmax; i++)
      p[i] = log(p[i]);

  return p;
}


// [[Rcpp::export]]
NumericVector cpp_qdweibull(NumericVector p,
                            NumericVector q, NumericVector beta,
                            bool lower_tail = true, bool log_prob = false) {

  if (is_true(any(q <= 0)) || is_true(any(q >= 1)))
    throw Rcpp::exception("Values of q should be 0 < q < 1.");
  if (is_true(any(beta <= 0)))
    throw Rcpp::exception("Values of beta should be > 0.");

  int n  = p.length();
  int nq = q.length();
  int nb = beta.length();
  int Nmax = Rcpp::max(IntegerVector::create(n, nq, nb));
  NumericVector x(Nmax);

  if (log_prob)
    for (int i = 0; i < n; i++)
      p[i] = exp(p[i]);

  if (!lower_tail)
    for (int i = 0; i < n; i++)
      p[i] = 1-p[i];

  for (int i = 0; i < Nmax; i++)
    x[i] = invcdf_dweibull(p[i % n], q[i % nq], beta[i % nb]);

  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rdweibull(int n,
                            NumericVector q, NumericVector beta) {

  if (is_true(any(q <= 0)) || is_true(any(q >= 1)))
    throw Rcpp::exception("Values of q should be 0 < q < 1.");
  if (is_true(any(beta <= 0)))
    throw Rcpp::exception("Values of beta should be > 0.");

  double u;
  int nq = q.length();
  int nb = beta.length();
  NumericVector x(n);

  for (int i = 0; i < n; i++) {
    u = R::runif(0, 1);
    x[i] = invcdf_dweibull(u, q[i % nq], beta[i % nb]);
  }

  return x;
}

