#ifndef EDCPP_SHARED_H
#define EDCPP_SHARED_H

// Basic functions

bool tol_equal(double x, double y);
bool isInteger(double x);

// Dealing with Inf

bool anyFinite(Rcpp::NumericVector x);
double finite_max(Rcpp::NumericVector x);

// Standard normal

double phi(double x);
double Phi(double x);
double InvPhi(double x);

// Factorial

double lfactorial(double x);
double factorial(double x);

// Random generation for Bernoulli

double rng_bernoulli(double p);

#endif
