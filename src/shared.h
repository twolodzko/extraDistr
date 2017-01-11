
#ifndef MIN_DIFF_EPS
#define MIN_DIFF_EPS 1e-8
#endif

#ifndef EDCPP_SHARED_H
#define EDCPP_SHARED_H

// Basic functions

bool isInteger(double x, bool warn = true);
double finite_max(Rcpp::NumericVector x);

// Random generation

double rng_unif();         // standard uniform
double rng_sign();         // Rademacher distribution


// ====================================================================
//                      Inline functions
// ====================================================================

inline bool tol_equal(double x, double y);
inline double phi(double x);
inline double Phi(double x);
inline double InvPhi(double x);
inline double factorial(double x);
inline double lfactorial(double x);


inline bool tol_equal(double x, double y) {
  return std::abs(x - y) < MIN_DIFF_EPS;
}

inline double phi(double x) {
  return R::dnorm(x, 0.0, 1.0, false);
}

inline double Phi(double x) {
  return R::pnorm(x, 0.0, 1.0, true, false);
}

inline double InvPhi(double x) {
  return R::qnorm(x, 0.0, 1.0, true, false);
}

inline double factorial(double x) {
  return R::gammafn(x + 1.0);
}

inline double lfactorial(double x) {
  return R::lgammafn(x + 1.0);
}


#endif
