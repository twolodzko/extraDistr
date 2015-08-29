#include <Rcpp.h>
using namespace Rcpp;


/*
*  Poisson-binomial distribution
*
*  Values:
*  x
*
*  Parameters:
*  0 < prob < 1
*
*/

int rng_bern(double p) {
  double u = R::runif(0, 1);
  if (u <= p)
    return 1;
  else
    return 0;
}

int rng_pbinom(NumericVector prob) {
  int k = prob.length();
  int x = 0;
  for (int i = 0; i < k; i++)
    x += rng_bern(prob[i]);
  return x;
}


// [[Rcpp::export]]
NumericVector cpp_rpbinom(int n,
                          NumericVector prob) {

  if (is_true(any(prob <= 0)) || is_true(any(prob >= 1)))
    throw Rcpp::exception("Values of theta should fit 0 < prob < 1.");

  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_pbinom(prob);

  return x;
}


// [[Rcpp::export]]
NumericVector sim_dpbinom(NumericVector prob,
                          int nsim = 10000) {

  int k = prob.length();
  NumericVector p_tab(k+1);

  for (int j = 0; j < nsim; j++) {
    double xsim = rng_pbinom(prob);
    p_tab[xsim] += 1.0/nsim;
  }

  return p_tab;
}


// [[Rcpp::export]]
NumericVector rna_ppbinom(NumericVector prob) {

  int k = prob.length();
  NumericVector p_tab(k+1);
  double mu = 0;
  double sigma = 0;
  double gamma = 0;

  for (int i = 0; i < k+1; i++) {
    double p_tmp = prob[i] * (1-prob[i]);
    mu += prob[i];
    sigma += p_tmp;
    gamma += p_tmp * (1-2*prob[i]);
  }
  sigma = sqrt(sigma);
  gamma *= pow(sigma, -3);

  for (int i = 0; i < k+1; i++) {
    double z = (i+0.5-mu)/sigma;
    double p_tmp = R::pnorm(z, 0, 1, true, false) + gamma*(1-pow(z, 2))
                   * R::dnorm(z, 0, 1, false)/6;
    if (p_tmp >= 0)
      p_tab[i] = p_tmp;
    else
      p_tab[i] = 0;
  }

  return p_tab;
}


/*
// [[Rcpp::export]]
NumericVector ra_dpbinom(NumericVector prob) {

  int k = prob.length();
  NumericVector p_tab(k+1), w(k), tj(k, 0.0);

  for (int i = 0; i < k; i++)
    w[i] = prob[i]/(1-prob[i]);

  for (int i = 0; i < k; i++)
    for (int j = 0; j < k; j++)
      tj[i] += pow(w[j], i+1);

  p_tab[0] = 1;
  for (int i = 0; i < k; i++)
    p_tab[0] *= (1-prob[i]);

  for (int i = 1; i < k+1; i++) {
    double p_tmp = 0.0;
    for (int j = 1; j < i+1; j++) {
      p_tmp += pow(-1, j-1) * p_tab[i-j] * tj[j-1];
    }
    p_tmp = p_tmp/i;

    if (p_tmp >= 0)
      p_tab[i] = p_tmp;
    else
      p_tab[i] = 0;
  }

  return p_tab;
}
*/

