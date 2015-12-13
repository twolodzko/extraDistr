#ifndef SHARED_H
#define SHARED_H

// Standard normal

double phi(double x);
double Phi(double x);
double InvPhi(double x);

// Error function

double erf(double x);
double erfc(double x);
double inv_erf(double x);

// Factorial

double lfactorial(double x);
double factorial(double x);

// Random generation for Bernoulli

int rng_bernoulli(double p);

#endif
