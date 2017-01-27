
#ifndef EDCPP_SHARED_H
#define EDCPP_SHARED_H

// Constants

#ifndef SQRT_2_PI
#define SQRT_2_PI 2.506628274631000241612
#endif

#ifndef PHI_0
#define PHI_0 0.3989422804014327028632
#endif

#ifndef MIN_DIFF_EPS
#define MIN_DIFF_EPS 1e-8
#endif

// MACROS

#define GETV(x, i)      x[i % x.length()]    // wrapped indexing of vector
#define GETM(x, i, j)   x(i % x.nrow(), j)   // wrapped indexing of matrix
#define VALID_PROB(p)   ((p >= 0.0) && (p <= 1.0))

// Basic functions

bool isInteger(double x, bool warn = true);
double finite_max(const Rcpp::NumericVector& x);
double rng_unif();         // standard uniform

// ====================================================================
//                      Inline functions
// ====================================================================

inline bool tol_equal(double x, double y);
inline double phi(double x);
inline double Phi(double x);
inline double InvPhi(double x);
inline double factorial(double x);
inline double lfactorial(double x);
inline double rng_sign();
inline double to_dbl(long int x);
inline long int to_int(double x);


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

inline double rng_sign() {
  double u = rng_unif();
  return (u > 0.5) ? 1.0 : -1.0;
}

inline double to_dbl(long int x) {
  return static_cast<double>(x);
}

inline long int to_int(double x) {
  if (R_FINITE(x) && x > std::numeric_limits<long int>::max())
    Rcpp::stop("reached largest integer which can be represented as <long int>");
  return static_cast<long int>(x);
}


#endif
