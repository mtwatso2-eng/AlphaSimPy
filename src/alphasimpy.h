#ifndef ALPHASIMPY_H
#define ALPHASIMPY_H

#include <armadillo>  // Use Armadillo directly instead of RcppArmadillo
#include <bitset>
// Don't include optimize.h for Python bindings - it uses Rcpp types
// #include "optimize.h"
#include "misc.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#endif
